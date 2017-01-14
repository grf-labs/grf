#include "ForestPredictor.h"
#include "TreePredictor.h"
#include "utility.h"

ForestPredictor::ForestPredictor(uint num_threads,
                                 std::shared_ptr<PredictionStrategy> prediction_strategy) :
    prediction_strategy(prediction_strategy) {
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

std::vector<std::vector<double>> ForestPredictor::predict(const Forest& forest, Data* prediction_data) {
  std::map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree =
      determine_terminal_node_IDs(forest, prediction_data, false);

  std::vector<std::vector<double>> predictions;
  predictions.reserve(prediction_data->getNumRows());

  size_t prediction_length = prediction_strategy->prediction_length();
  for (size_t sampleID = 0; sampleID < prediction_data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;

    // Calculate the weights of neighboring samples.
    for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
      std::vector<size_t> terminal_node_IDs = terminal_node_IDs_by_tree.at(tree_index);
      add_sample_weights(sampleID, tree, weights_by_sampleID, terminal_node_IDs);
    }
    normalize_sample_weights(weights_by_sampleID);

    std::vector<double> prediction = prediction_strategy->predict(weights_by_sampleID,
                                                                  forest.get_observations());
    if (prediction.size() != prediction_length) {
      throw std::runtime_error("Prediction for sample " + std::to_string(sampleID) +
                               " did not have the expected length.");
    }
    predictions.push_back(prediction);
  }

  return predictions;
}

std::vector<std::vector<double>> ForestPredictor::predict_oob(const Forest& forest, Data* original_data) {
  std::map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree =
      determine_terminal_node_IDs(forest, original_data, true);

  std::vector<std::vector<double>> predictions;
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
    std::unordered_map<size_t, double> weights_by_sampleID;

    // If this sample was never out-of-bag, then return placeholder predictions.
    auto tree_range = trees_by_oob_samples.equal_range(sampleID);
    if (tree_range.first == tree_range.second) {
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(temp);
      continue;
    }

    // Calculate the weights of neighboring samples.
    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      size_t tree_index = it->second;
      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
      add_sample_weights(sampleID, tree, weights_by_sampleID, terminal_node_IDs_by_tree.at(tree_index));
    }

    // If this sample has no neighbors, then return placeholder predictions. Note
    // that this can only occur when honesty is enabled, and is expected to be rare.
    if (weights_by_sampleID.empty()) {
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(temp);
      continue;
    }

    normalize_sample_weights(weights_by_sampleID);

    std::vector<double> prediction = prediction_strategy->predict(weights_by_sampleID,
                                                                  forest.get_observations());
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
                                           const Data *prediction_data,
                                           bool oob_prediction,
                                           std::promise<std::map<size_t, std::vector<size_t>>> promise) {
  std::map<size_t, std::vector<size_t>> terminal_nodeIDs_by_tree;
  TreePredictor tree_predictor;

  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      std::shared_ptr<Tree> tree = forest.get_trees()[i];

      std::vector<size_t> sampleIDs;
      if (oob_prediction) {
        sampleIDs = tree->get_oob_sampleIDs();
      }

      std::vector<size_t> terminal_nodeIDs = tree_predictor.get_terminal_nodeIDs(tree,
          prediction_data,
          sampleIDs);
      terminal_nodeIDs_by_tree[i] = terminal_nodeIDs;
    }
  }
  promise.set_value(terminal_nodeIDs_by_tree);
}

void ForestPredictor::add_sample_weights(size_t test_sampleID,
                                         std::shared_ptr<Tree> tree,
                                         std::unordered_map<size_t, double> &weights_by_sampleID,
                                         std::vector<size_t> terminal_node_IDs) {
  size_t nodeID = terminal_node_IDs[test_sampleID];
  std::vector<size_t> sampleIDs = tree->get_terminal_nodeIDs()[nodeID];
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestPredictor::normalize_sample_weights(std::unordered_map<size_t, double> &weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}
