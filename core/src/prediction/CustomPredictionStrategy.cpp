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

#include "CustomPredictionStrategy.h"

const std::size_t CustomPredictionStrategy::OUTCOME = 0;

size_t CustomPredictionStrategy::prediction_length() {
  return 1;
}

Prediction CustomPredictionStrategy::predict(size_t sampleID,
                                             const std::vector<double>& averages,
                                             const std::unordered_map<size_t, double>& weights_by_sampleID,
                                             const Observations& observations) {
  std::vector<double> prediction = { 0.0 };
  return Prediction(prediction);
}

Prediction CustomPredictionStrategy::predict_with_variance(
    size_t sampleID,
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations,
    uint ci_group_size) {
  throw std::runtime_error("Variance estimates are not yet implemented.");
}

bool CustomPredictionStrategy::requires_leaf_sampleIDs() {
  return true;
}

PredictionValues CustomPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations) {
  return PredictionValues();
}
