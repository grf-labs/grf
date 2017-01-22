#include "ForestPredictor.h"
#include "Prediction.h"
#include "TreePredictor.h"
#include "utility.h"

ForestPredictor::ForestPredictor(uint num_threads,
                                 uint ci_group_size,
                                 std::shared_ptr<PredictionStrategy> prediction_strategy) :
    ci_group_size(ci_group_size), prediction_strategy(prediction_strategy) {
  if (num_threads == DEFAULT_NUM_THREADS) {
    this->num_threads = std::thread::hardware_concurrency();
  } else {
    this->num_threads = num_threads;
  }
}

std::map<size_t, std::vector<size_t>> ForestPredictor::determine_terminal_node_IDs(
    const Forest& forest,
    Data* data,
    bool oob_prediction) {
  std::map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree;

  size_t num_trees = forest.get_trees().size();
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::vector<std::future<
      std::map<size_t, std::vector<size_t>>>> futures;

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::map<size_t, std::vector<size_t>>> promise;
    std::future<std::map<size_t, std::vector<size_t>>> future = promise.get_future();
    threads.push_back(std::thread(&ForestPredictor::predictTreesInThread,
                                  this,
                                  i,
                                  forest,
                                  data,
                                  oob_prediction,
                                  std::move(promise)));
    futures.push_back(std::move(future));
  }

  for (auto &future : futures) {
    future.wait();
    std::map<size_t, std::vector<size_t>> terminal_nodeIDs = future.get();
    terminal_node_IDs_by_tree.insert(terminal_nodeIDs.begin(), terminal_nodeIDs.end());
  }

  for (auto& thread : threads) {
    thread.join();
  }
  return terminal_node_IDs_by_tree;
};

std::vector<Prediction> ForestPredictor::predict(const Forest& forest, Data* prediction_data) {
  std::map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree =
      determine_terminal_node_IDs(forest, prediction_data, false);

  std::vector<Prediction> predictions;
  predictions.reserve(prediction_data->getNumRows());

  size_t prediction_length = prediction_strategy->prediction_length();
  for (size_t sampleID = 0; sampleID < prediction_data->getNumRows(); ++sampleID) {
    std::map<std::string, double> average_prediction_values;
    std::unordered_map<size_t, double> weights_by_sampleID;
    std::vector<std::vector<size_t>> leaf_sampleIDs;

    // Average the precomputed prediction values across the leaf nodes this sample falls
    // in, and create a list of weighted neighbors if the prediction strategy requires it.
    uint num_nonempty_leaves = 0;
    for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];

      std::vector<size_t> terminal_node_IDs = terminal_node_IDs_by_tree.at(tree_index);
      size_t nodeID = terminal_node_IDs.at(sampleID);
      const std::vector<size_t>& sampleIDs = tree->get_leaf_nodeIDs()[nodeID];

      // hackhackhack
      if (ci_group_size == 1) {
        if (sampleIDs.empty()) {
          continue;
        }

        num_nonempty_leaves++;
        add_prediction_values(nodeID, tree->get_prediction_values(), average_prediction_values);
        if (prediction_strategy->requires_leaf_sampleIDs()) {
          add_sample_weights(sampleIDs, weights_by_sampleID);
        }
      } else {
        leaf_sampleIDs.push_back(sampleIDs);
      }
    }

    if (ci_group_size == 1) {
      if (num_nonempty_leaves == 0) {
        // If this sample has no neighbors, then return placeholder predictions. Note
        // that this can only occur when honesty is enabled, and is expected to be rare.
        std::vector<double> temp(prediction_length, NAN);
        predictions.push_back(temp);
        continue;
      } else {
        normalize_prediction_values(num_nonempty_leaves, average_prediction_values);
        normalize_sample_weights(weights_by_sampleID);
      }
    }

    Prediction prediction = ci_group_size == 1 ?
        prediction_strategy->predict(average_prediction_values, weights_by_sampleID, forest.get_observations()) :
        prediction_strategy->predict_with_variance(leaf_sampleIDs, forest.get_observations(), ci_group_size);

    if (prediction.size() != prediction_length) {
      throw std::runtime_error("Prediction for sample " + std::to_string(sampleID) +
                               " did not have the expected length.");
    }

    predictions.push_back(prediction);
  }

  return predictions;
}

std::vector<Prediction> ForestPredictor::predict_oob(const Forest& forest, Data* original_data) {
  std::map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree =
      determine_terminal_node_IDs(forest, original_data, true);

  std::vector<Prediction> predictions;
  predictions.reserve(original_data->getNumRows());

  size_t num_trees = forest.get_trees().size();
  std::unordered_multimap<size_t, size_t> trees_by_oob_samples;
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sampleID : forest.get_trees()[tree_idx]->get_oob_sampleIDs()) {
      trees_by_oob_samples.insert(std::pair<size_t, size_t>(sampleID, tree_idx));
    }
  }

  size_t prediction_length = prediction_strategy->prediction_length();
  for (size_t sampleID = 0; sampleID < original_data->getNumRows(); ++sampleID) {
    std::map<std::string, double> average_prediction_values;
    std::unordered_map<size_t, double> weights_by_sampleID;

    // If this sample was never out-of-bag, then return placeholder predictions.
    auto tree_range = trees_by_oob_samples.equal_range(sampleID);
    if (tree_range.first == tree_range.second) {
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(temp);
      continue;
    }

    // Average the precomputed prediction values across the leaf nodes this sample falls
    // in, and create a list of weighted neighbors if the prediction strategy requires it.
    uint num_nonempty_leaves = 0;
    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      size_t tree_index = it->second;
      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];

      std::vector<size_t> terminal_node_IDs = terminal_node_IDs_by_tree.at(tree_index);
      size_t nodeID = terminal_node_IDs.at(sampleID);
      const std::vector<size_t>& sampleIDs = tree->get_leaf_nodeIDs()[nodeID];

      if (sampleIDs.empty()) {
        continue;
      }

      num_nonempty_leaves++;
      add_prediction_values(nodeID, tree->get_prediction_values(), average_prediction_values);
      if (prediction_strategy->requires_leaf_sampleIDs()) {
        add_sample_weights(sampleIDs, weights_by_sampleID);
      }
    }

    if (ci_group_size == 1) {
      if (num_nonempty_leaves == 0) {
        // If this sample has no neighbors, then return placeholder predictions. Note
        // that this can only occur when honesty is enabled, and is expected to be rare.
        std::vector<double> temp(prediction_length, NAN);
        predictions.push_back(temp);
        continue;
      } else {
        normalize_prediction_values(num_nonempty_leaves, average_prediction_values);
        normalize_sample_weights(weights_by_sampleID);
      }
    }

    Prediction prediction = prediction_strategy->predict(average_prediction_values,
        weights_by_sampleID, forest.get_observations());

    if (prediction.size() != prediction_length) {
      throw std::runtime_error("Prediction for sample " + std::to_string(sampleID) +
          " did not have the expected length.");
    }
    predictions.push_back(prediction);
  }

  return predictions;
}

void ForestPredictor::predictTreesInThread(uint thread_idx,
                                           const Forest& forest,
                                           Data* prediction_data,
                                           bool oob_prediction,
                                           std::promise<std::map<size_t, std::vector<size_t>>> promise) {
  std::map<size_t, std::vector<size_t>> terminal_nodeIDs_by_tree;
  TreePredictor tree_predictor;

  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      std::shared_ptr<Tree> tree = forest.get_trees()[i];

      const std::vector<size_t>& sampleIDs = oob_prediction ? tree->get_oob_sampleIDs() : std::vector<size_t>();

      std::vector<size_t> terminal_nodeIDs = tree_predictor.get_terminal_nodeIDs(tree,
          prediction_data,
          sampleIDs);
      terminal_nodeIDs_by_tree[i] = terminal_nodeIDs;
    }
  }
  promise.set_value(terminal_nodeIDs_by_tree);
}

void ForestPredictor::add_prediction_values(size_t nodeID,
                                            const PredictionValues &prediction_values,
                                            std::map<std::string, double>& average_prediction_values) {
  auto& values_by_type = prediction_values.get_values_by_type();
  for (auto it = values_by_type.begin(); it != values_by_type.end(); it++) {
    std::string type = it->first;
    average_prediction_values[type] += std::isnan(it->second[nodeID]) ? 0 : it->second[nodeID];
  }
}

void ForestPredictor::add_sample_weights(const std::vector<size_t>& sampleIDs,
                                         std::unordered_map<size_t, double>& weights_by_sampleID) {
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestPredictor::normalize_prediction_values(size_t num_leaf_nodes,
                                                  std::map<std::string, double>& average_prediction_values) {
  for (auto it = average_prediction_values.begin(); it != average_prediction_values.end(); ++it) {
    it->second /= num_leaf_nodes;
  }
}

void ForestPredictor::normalize_sample_weights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}
