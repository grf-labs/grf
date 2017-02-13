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

const std::string RegressionPredictionStrategy::OUTCOME = "outcome";

size_t RegressionPredictionStrategy::prediction_length() {
    return 1;
}

Prediction RegressionPredictionStrategy::predict(const std::map<std::string, double>& averages,
                                                 const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                 const Observations& observations) {
  std::vector<double> predictions = { averages.at(OUTCOME) };
  return Prediction(predictions);
}

Prediction RegressionPredictionStrategy::predict_with_variance(
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

  std::map<std::string, std::vector<double>> values;
  std::vector<double>& averages = values[OUTCOME];
  averages.resize(num_leaves);

  for (int i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_sampleIDs.at(i);

    if (leaf_node.empty()) {
      averages[i] = NAN;
      continue;
    }

    double average = 0.0;
    for (auto& sampleID : leaf_node) {
      average += observations.get(Observations::OUTCOME).at(sampleID);
    }
    averages[i] = average / leaf_node.size();
  }

  return PredictionValues(values, num_leaves);
}

