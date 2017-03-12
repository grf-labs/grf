/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "ForestPredictor.h"
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

std::vector<std::vector<size_t>> ForestPredictor::determine_terminal_node_IDs(
    const Forest& forest,
    Data* data,
    bool oob_prediction) {
  size_t num_trees = forest.get_trees().size();

  std::vector<std::vector<size_t>> all_terminal_node_IDs;
  all_terminal_node_IDs.reserve(num_trees);

  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, num_trees - 1, num_threads);

  std::vector<std::future<
      std::vector<std::vector<size_t>>>> futures;
  futures.reserve(num_threads);

  for (uint i = 0; i < num_threads; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_trees_batch = thread_ranges[i + 1] - start_index;
    futures.push_back(std::async(std::launch::async,
                                 &ForestPredictor::predict_batch,
                                 this,
                                 start_index,
                                 num_trees_batch,
                                 forest,
                                 data,
                                 oob_prediction));
  }

  for (auto& future : futures) {
    std::vector<std::vector<size_t>> terminal_node_IDs = future.get();
    all_terminal_node_IDs.insert(all_terminal_node_IDs.end(),
                                 terminal_node_IDs.begin(),
                                 terminal_node_IDs.end());
  }

  return all_terminal_node_IDs;
};

std::vector<Prediction> ForestPredictor::standard_predict(const Forest& forest, Data* prediction_data) {
  size_t num_trees = forest.get_trees().size();
  size_t num_samples = prediction_data->get_num_rows();

  std::vector<std::vector<size_t>> terminal_nodeIDs_by_tree =
      determine_terminal_node_IDs(forest, prediction_data, false);

  // Collect the average prediction values across trees.
  std::vector<std::vector<double>> average_prediction_values(num_samples);
  std::vector<size_t> num_nonempty_leaves(num_samples);

  for (size_t tree_index = 0; tree_index < num_trees; ++tree_index) {
    std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
    const PredictionValues &prediction_values = tree->get_prediction_values();

    const std::vector<size_t>& terminal_node_IDs = terminal_nodeIDs_by_tree.at(tree_index);

    for (size_t sampleID = 0; sampleID < prediction_data->get_num_rows(); sampleID++) {
      size_t nodeID = terminal_node_IDs.at(sampleID);

      if (!prediction_values.empty(nodeID)) {
        add_prediction_values(nodeID, prediction_values, average_prediction_values[sampleID]);
        num_nonempty_leaves[sampleID]++;
      }
    }
  }

  // Normalize the prediction values, and make a prediction for each sample.
  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);
  size_t prediction_length = prediction_strategy->prediction_length();

  for (size_t sampleID = 0; sampleID < prediction_data->get_num_rows(); sampleID++) {
    size_t num_leaves = num_nonempty_leaves[sampleID];
    if (num_leaves == 0) {
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(Prediction(temp));
    } else {
      std::vector<double> &average_prediction_value = average_prediction_values[sampleID];
      normalize_prediction_values(num_leaves, average_prediction_value);
      predictions.push_back(prediction_strategy->predict(average_prediction_value,
                                                         std::unordered_map<size_t, double>(),
                                                         forest.get_observations()));
    }
  }
  return predictions;
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest, Data* prediction_data) {
  if (!prediction_strategy->requires_leaf_sampleIDs() && ci_group_size == 1) {
    return standard_predict(forest, prediction_data);
  }

  size_t num_trees = forest.get_trees().size();
  std::vector<std::vector<size_t>> all_terminal_node_IDs =
      determine_terminal_node_IDs(forest, prediction_data, false);

  std::vector<Prediction> predictions;
  predictions.reserve(prediction_data->get_num_rows());

  size_t prediction_length = prediction_strategy->prediction_length();

  for (size_t sampleID = 0; sampleID < prediction_data->get_num_rows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;

    std::vector<std::vector<size_t>> leaf_sampleIDs;
    if (ci_group_size > 1) {
      leaf_sampleIDs.reserve(num_trees);
    }

    // Average the precomputed prediction values across the leaf nodes this sample falls
    // in, and create a list of weighted neighbors if the prediction strategy requires it.
    uint num_nonempty_leaves = 0;
    for (size_t tree_index = 0; tree_index < num_trees; ++tree_index) {
      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];

      const std::vector<size_t>& terminal_node_IDs = all_terminal_node_IDs.at(tree_index);
      size_t nodeID = terminal_node_IDs.at(sampleID);
      const std::vector<size_t>& sampleIDs = tree->get_leaf_nodeIDs()[nodeID];

      if (!sampleIDs.empty()) {
        num_nonempty_leaves++;
      }

      // hackhackhack
      if (ci_group_size == 1) {
        if (!sampleIDs.empty()) {
          add_sample_weights(sampleIDs, weights_by_sampleID);
        }
      } else {
        leaf_sampleIDs.push_back(sampleIDs);
      }
    }

    if (num_nonempty_leaves == 0) {
      // If this sample has no neighbors, then return placeholder predictions. Note
      // that this can only occur when honesty is enabled, and is expected to be rare.
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(temp);
      continue;
    } else {
      if (ci_group_size == 1) {
        normalize_sample_weights(weights_by_sampleID);
      }
    }

    Prediction prediction = ci_group_size == 1 ?
        prediction_strategy->predict(std::vector<double>(), weights_by_sampleID, forest.get_observations()) :
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
  std::vector<std::vector<size_t>> terminal_node_IDs_by_tree =
      determine_terminal_node_IDs(forest, original_data, true);

  std::vector<Prediction> predictions;
  predictions.reserve(original_data->get_num_rows());

  size_t num_trees = forest.get_trees().size();
  bool trees_by_oob_samples[original_data->get_num_rows()][num_trees];
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sampleID : forest.get_trees()[tree_idx]->get_oob_sampleIDs()) {
      trees_by_oob_samples[sampleID][tree_idx] = true;
    }
  }

  size_t prediction_length = prediction_strategy->prediction_length();

  for (size_t sampleID = 0; sampleID < original_data->get_num_rows(); ++sampleID) {
    std::vector<double> average_prediction_values;
    std::unordered_map<size_t, double> weights_by_sampleID;

    // Average the precomputed prediction values across the leaf nodes this sample falls
    // in, and create a list of weighted neighbors if the prediction strategy requires it.
    uint num_nonempty_leaves = 0;
    for (size_t t = 0; t < num_trees; t++) {
      if (!trees_by_oob_samples[sampleID][t]) {
        continue;
      }

      std::shared_ptr<Tree> tree = forest.get_trees()[t];
      const std::vector<size_t>& terminal_node_IDs = terminal_node_IDs_by_tree.at(t);
      size_t nodeID = terminal_node_IDs.at(sampleID);
      const std::vector<size_t>& sampleIDs = tree->get_leaf_nodeIDs()[nodeID];

      if (!sampleIDs.empty()) {
        num_nonempty_leaves++;
        add_prediction_values(nodeID, tree->get_prediction_values(), average_prediction_values);
        if (prediction_strategy->requires_leaf_sampleIDs()) {
          add_sample_weights(sampleIDs, weights_by_sampleID);
        }
      }
    }

    if (num_nonempty_leaves == 0) {
      // If this sample has no neighbors, then return placeholder predictions. This could
      // occur because the sample is never OOB, or in rare cases when honesty is enabled.
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(temp);
      continue;
    } else {
      normalize_prediction_values(num_nonempty_leaves, average_prediction_values);
      normalize_sample_weights(weights_by_sampleID);
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

std::vector<std::vector<size_t>> ForestPredictor::predict_batch(
    size_t start,
    size_t num_trees,
    const Forest& forest,
    Data* prediction_data,
    bool oob_prediction) {
  std::vector<std::vector<size_t>> all_terminal_node_IDs(num_trees);
  for (size_t i = 0; i < num_trees; i++) {
    std::shared_ptr<Tree> tree = forest.get_trees()[start + i];

    const std::vector<size_t> &sampleIDs = oob_prediction ? tree->get_oob_sampleIDs() : std::vector<size_t>();
    std::vector<size_t> terminal_node_IDs = tree->find_leaf_nodeIDs(prediction_data, sampleIDs);
    all_terminal_node_IDs[i] = terminal_node_IDs;
  }

  return all_terminal_node_IDs;
}

void ForestPredictor::add_prediction_values(size_t nodeID,
                                            const PredictionValues& prediction_values,
                                            std::vector<double>& average_prediction_values) {
  if (average_prediction_values.empty()) {
    average_prediction_values.resize(prediction_values.get_num_types());
  }

  for (size_t type = 0; type < prediction_values.get_num_types(); type++) {
    average_prediction_values[type] += prediction_values.get(nodeID, type);
  }
}

void ForestPredictor::add_sample_weights(const std::vector<size_t>& sampleIDs,
                                         std::unordered_map<size_t, double>& weights_by_sampleID) {
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto& sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestPredictor::normalize_prediction_values(size_t num_leaf_nodes,
                                                  std::vector<double>& average_prediction_values) {
  for (double& value : average_prediction_values) {
    value /= num_leaf_nodes;
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
