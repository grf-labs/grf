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

#include "prediction/collector/DefaultPredictionCollector.h"

DefaultPredictionCollector::DefaultPredictionCollector(std::shared_ptr<DefaultPredictionStrategy> prediction_strategy):
    prediction_strategy(prediction_strategy) {}

std::vector<Prediction> DefaultPredictionCollector::collect_predictions(
    const Forest& forest,
    Data* prediction_data,
    std::vector<std::vector<size_t>> leaf_nodes_by_tree,
    bool oob_prediction) {

  size_t num_samples = prediction_data->get_num_rows();
  size_t num_trees = forest.get_trees().size();

  std::vector<std::vector<bool>> trees_by_samples(num_samples, std::vector<bool>(num_trees));
  if (oob_prediction) {
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      for (size_t sampleID : forest.get_trees()[tree_idx]->get_oob_sampleIDs()) {
        trees_by_samples[sampleID][tree_idx] = true;
      }
    }
  }

  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);
  size_t prediction_length = prediction_strategy->prediction_length();

  for (size_t sampleID = 0; sampleID < num_samples; ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;


    // Create a list of weighted neighbors for this sample.
    uint num_nonempty_leaves = 0;
    for (size_t tree_index = 0; tree_index < num_trees; ++tree_index) {
      if (oob_prediction && !trees_by_samples[sampleID][tree_index]) {
        continue;
      }

      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
      const std::vector<size_t> &leaf_node_IDs = leaf_nodes_by_tree.at(tree_index);

      size_t nodeID = leaf_node_IDs.at(sampleID);
      const std::vector<size_t> &sampleIDs = tree->get_leaf_nodeIDs()[nodeID];
      if (!sampleIDs.empty()) {
        num_nonempty_leaves++;
        add_sample_weights(sampleIDs, weights_by_sampleID);
      }
    }

    if (num_nonempty_leaves == 0) {
      // If this sample has no neighbors, then return placeholder predictions. Note
      // that this can only occur when honesty is enabled, and is expected to be rare.
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(temp);
      continue;
    }

    normalize_sample_weights(weights_by_sampleID);

    Prediction prediction = prediction_strategy->predict(
        sampleID, weights_by_sampleID, forest.get_observations());
    if (prediction.size() != prediction_length) {
      throw std::runtime_error("Prediction for sample " + std::to_string(sampleID) +
                               " did not have the expected length.");
    }
    predictions.push_back(prediction);
  }
  return predictions;
}

void DefaultPredictionCollector::add_sample_weights(const std::vector<size_t>& sampleIDs,
                                          std::unordered_map<size_t, double>& weights_by_sampleID) {
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto& sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void DefaultPredictionCollector::normalize_sample_weights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}
