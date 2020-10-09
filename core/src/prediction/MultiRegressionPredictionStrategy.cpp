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

#include "prediction/MultiRegressionPredictionStrategy.h"

#include "Eigen/Dense"

namespace grf {

MultiRegressionPredictionStrategy::MultiRegressionPredictionStrategy(size_t num_outcomes) {
  this->num_outcomes = num_outcomes;
}

size_t MultiRegressionPredictionStrategy::prediction_length() const {
  return num_outcomes;
}

std::vector<double> MultiRegressionPredictionStrategy::predict(const std::vector<double>& average) const {
  return average;
}

std::vector<double> MultiRegressionPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {
  return { 0.0 };
}


size_t MultiRegressionPredictionStrategy::prediction_value_length() const {
  return 1;
}

PredictionValues MultiRegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_samples.at(i);
    if (leaf_node.empty()) {
      continue;
    }

    Eigen::VectorXd sum = Eigen::VectorXd::Zero(num_outcomes);
    double weight = 0.0;
    for (auto& sample : leaf_node) {
      sum += data.get_weight(sample) * data.get_outcomes(sample);
      weight += data.get_weight(sample);
    }

    // if total weight is very small, treat the leaf as empty
    if (std::abs(weight) <= 1e-16) {
      continue;
    }
    sum /= weight;

    values[i] = std::vector<double> (sum.data(), sum.data() + num_outcomes);
  }

  return PredictionValues(values, num_outcomes);
}

std::vector<std::pair<double, double>> MultiRegressionPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
