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

OptimizedPredictionCollector::OptimizedPredictionCollector(std::shared_ptr<OptimizedPredictionStrategy> strategy,
                                                           uint ci_group_size):
    strategy(strategy),
    ci_group_size(ci_group_size) {}

std::vector<Prediction> OptimizedPredictionCollector::collect_predictions(const Forest &forest,
                                                                          Data *prediction_data,
                                                                          std::vector<std::vector<size_t>> leaf_nodes_by_tree,
                                                                          std::vector<std::vector<bool>> trees_by_sample) {
  size_t num_trees = forest.get_trees().size();
  size_t num_samples = prediction_data->get_num_rows();

  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  for (size_t sampleID = 0; sampleID < num_samples; ++sampleID) {
    std::vector<double> average_value;
    std::vector<std::vector<double>> leaf_values;
    if (ci_group_size > 1) {
      leaf_values.resize(num_trees);
    }

    // Create a list of weighted neighbors for this sample.
    uint num_leaves = 0;
    for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
      if (!trees_by_sample.empty() && !trees_by_sample[sampleID][tree_index]) {
        continue;
      }

      const std::vector<size_t>& leaf_node_IDs = leaf_nodes_by_tree.at(tree_index);
      size_t nodeID = leaf_node_IDs.at(sampleID);

      std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
      const PredictionValues& prediction_values = tree->get_prediction_values();

      if (!prediction_values.empty(nodeID)) {
        num_leaves++;
        add_prediction_values(nodeID, prediction_values, average_value);
        if (ci_group_size > 1) {
          leaf_values[tree_index] = prediction_values.get_values(nodeID);
        }
      }
    }

    // If this sample has no neighbors, then return placeholder predictions. Note
    // that this can only occur when honesty is enabled, and is expected to be rare.
    if (num_leaves == 0) {
      std::vector<double> temp(strategy->prediction_length(), NAN);
      predictions.push_back(Prediction(temp));
      continue;
    }

    normalize_prediction_values(num_leaves, average_value);

    std::vector<double> point_prediction = strategy->predict(average_value);
    std::vector<double> variance_estimate;
    if (ci_group_size > 1) {
      variance_estimate = strategy->compute_variance(average_value, leaf_values, ci_group_size);
    }

    Prediction prediction(point_prediction, variance_estimate);
    validate_prediction(sampleID, prediction);
    predictions.push_back(prediction);
  }
  return predictions;
}

void OptimizedPredictionCollector::add_prediction_values(size_t nodeID,
    const PredictionValues& prediction_values,
    std::vector<double>& combined_average) {
  if (combined_average.empty()) {
    combined_average.resize(prediction_values.get_num_types());
  }

  for (size_t type = 0; type < prediction_values.get_num_types(); ++type) {
    combined_average[type] += prediction_values.get(nodeID, type);
  }
}

void OptimizedPredictionCollector::normalize_prediction_values(size_t num_leaves,
    std::vector<double>& combined_average) {
  for (double& value : combined_average) {
    value /= num_leaves;
  }
}

void OptimizedPredictionCollector::validate_prediction(size_t sampleID, Prediction prediction) {
  size_t prediction_length = strategy->prediction_length();
  if (prediction.size() != prediction_length) {
    throw std::runtime_error("Prediction for sample " + std::to_string(sampleID) +
                             " did not have the expected length.");
  }
}