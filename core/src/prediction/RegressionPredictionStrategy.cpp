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

#include <cmath>
#include <string>
#include "RegressionPredictionStrategy.h"

const size_t RegressionPredictionStrategy::OUTCOME = 0;

size_t RegressionPredictionStrategy::prediction_length() {
    return 1;
}

Prediction RegressionPredictionStrategy::predict(size_t sampleID,
                                                 const std::vector<double>& averages,
                                                 const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                 const Observations& observations) {
  std::vector<double> predictions = { averages.at(OUTCOME) };
  return Prediction(predictions);
}

Prediction RegressionPredictionStrategy::predict_with_variance(
    size_t sampleID,
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations,
    uint ci_group_size) {
  throw std::runtime_error("Variance estimates are not yet implemented.");
}

bool RegressionPredictionStrategy::requires_leaf_sampleIDs() {
  return false;
}

PredictionValues RegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations) {
  size_t num_leaves = leaf_sampleIDs.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_sampleIDs.at(i);
    if (leaf_node.empty()) {
      continue;
    }

    std::vector<double>& averages = values[i];
    averages.resize(1);

    double average = 0.0;
    for (auto& sampleID : leaf_node) {
      average += observations.get(Observations::OUTCOME, sampleID);
    }
    averages[OUTCOME] = average / leaf_node.size();
  }

  return PredictionValues(values, num_leaves, 1);
}

