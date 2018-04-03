/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "prediction/collector/DefaultPredictionCollector.h"

DefaultPredictionCollector::DefaultPredictionCollector(std::shared_ptr<DefaultPredictionStrategy> strategy):
    strategy(strategy) {}

std::vector<Prediction> DefaultPredictionCollector::collect_predictions(
    const Forest& forest,
    Data* prediction_data,
    const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
    const std::vector<std::vector<bool>>& valid_trees_by_sample,
    bool estimate_error) {

  size_t num_samples = prediction_data->get_num_rows();
  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  for (size_t sample = 0; sample < num_samples; sample++) {
    std::unordered_map<size_t, double> weights_by_sample = weight_computer.compute_weights(
        sample, forest, leaf_nodes_by_tree, valid_trees_by_sample);
    Prediction prediction = strategy->predict(sample, weights_by_sample, forest.get_observations());

    validate_prediction(sample, prediction);
    predictions.push_back(prediction);
  }
  return predictions;
}

void DefaultPredictionCollector::validate_prediction(size_t sample, Prediction prediction) {
  size_t prediction_length = strategy->prediction_length();
  if (prediction.size() != prediction_length) {
    throw std::runtime_error("Prediction for sample " + std::to_string(sample) +
                             " did not have the expected length.");
  }
}
