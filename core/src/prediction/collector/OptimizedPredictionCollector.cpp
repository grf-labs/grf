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

#include "prediction/collector/OptimizedPredictionCollector.h"

OptimizedPredictionCollector::OptimizedPredictionCollector(std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy):
    prediction_strategy(prediction_strategy) {}

std::vector<Prediction> OptimizedPredictionCollector::collect_predictions(const Forest &forest,
                                                                          Data *prediction_data,
                                                                          std::vector<std::vector<size_t>> leaf_nodes_by_tree,
                                                                          std::vector<std::vector<bool>> trees_by_sample) {
  size_t num_trees = forest.get_trees().size();
  size_t num_predictions = prediction_data->get_num_rows();
  bool oob_prediction = !trees_by_sample.empty();

  // Collect the average prediction values across trees.
  std::vector<std::vector<double>> average_prediction_values(num_predictions);
  std::vector<size_t> num_nonempty_leaves(num_predictions);

  for (size_t tree_index = 0; tree_index < num_trees; ++tree_index) {
    std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
    const PredictionValues &prediction_values = tree->get_prediction_values();
    const std::vector<size_t>& leaf_nodes = leaf_nodes_by_tree.at(tree_index);

    size_t num_samples = oob_prediction ? tree->get_oob_sampleIDs().size() : num_predictions;
    for (size_t i = 0; i < num_samples; ++i) {
      size_t sampleID = oob_prediction ? tree->get_oob_sampleIDs()[i] : i;
      size_t nodeID = leaf_nodes.at(sampleID);

      if (!prediction_values.empty(nodeID)) {
        add_prediction_values(nodeID, prediction_values, average_prediction_values[sampleID]);
        num_nonempty_leaves[sampleID]++;
      }
    }
  }

  // Normalize the prediction values, and make a prediction for each sample.
  std::vector<Prediction> predictions;
  predictions.reserve(num_predictions);
  size_t prediction_length = prediction_strategy->prediction_length();

  for (size_t sampleID = 0; sampleID < prediction_data->get_num_rows(); ++sampleID) {
    size_t num_leaves = num_nonempty_leaves[sampleID];

    if (num_leaves == 0) {
      // If this sample has no neighbors, then return placeholder predictions. Note
      // that this can only occur when honesty is enabled, and is expected to be rare.
      std::vector<double> temp(prediction_length, NAN);
      predictions.push_back(Prediction(temp));
    } else {
      std::vector<double> &average_prediction_value = average_prediction_values[sampleID];
      normalize_prediction_values(num_leaves, average_prediction_value);
      predictions.push_back(prediction_strategy->predict(average_prediction_value));
    }
  }
  return predictions;
}

void OptimizedPredictionCollector::add_prediction_values(size_t nodeID,
                                              const PredictionValues& prediction_values,
                                              std::vector<double>& average_prediction_values) {
  if (average_prediction_values.empty()) {
    average_prediction_values.resize(prediction_values.get_num_types());
  }

  for (size_t type = 0; type < prediction_values.get_num_types(); ++type) {
    average_prediction_values[type] += prediction_values.get(nodeID, type);
  }
}

void OptimizedPredictionCollector::normalize_prediction_values(size_t num_leaf_nodes,
                                                    std::vector<double>& average_prediction_values) {
  for (double& value : average_prediction_values) {
    value /= num_leaf_nodes;
  }
}