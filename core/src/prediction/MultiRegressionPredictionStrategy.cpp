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
  this->num_types = 1 + num_outcomes;
  this->weight_index = num_outcomes;
}

size_t MultiRegressionPredictionStrategy::prediction_length() const {
  return num_outcomes;
}

std::vector<double> MultiRegressionPredictionStrategy::predict(const std::vector<double>& average) const {
  std::vector<double> predictions;
  predictions.reserve(num_outcomes);
  double weight_bar = average[weight_index];
  for (size_t j = 0; j < num_outcomes; j++) {
    predictions.push_back(average[j] / weight_bar);
  }

  return predictions;
}

std::vector<double> MultiRegressionPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {
  return { 0.0 };
}

size_t MultiRegressionPredictionStrategy::prediction_value_length() const {
  return num_types;
}

PredictionValues MultiRegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_samples.at(i);
    size_t num_samples = leaf_node.size();
    if (num_samples == 0) {
      continue;
    }

    Eigen::VectorXd sum = Eigen::VectorXd::Zero(num_outcomes);
    double sum_weight = 0.0;
    for (auto& sample : leaf_node) {
      double weight = data.get_weight(sample);
      sum += weight * data.get_outcomes(sample);
      sum_weight += weight;
    }
    // if total weight is very small, treat the leaf as empty
    if (std::abs(sum_weight) <= 1e-16) {
      continue;
    }

    // store sufficient statistics in order
    // {outcome_1, ..., outcome_M, weight_sum}
    std::vector<double>& value = values[i];
    value.reserve(num_types);
    for (size_t j = 0; j < num_outcomes; j++) {
      value.push_back(sum[j] / num_samples);
    }
    value.push_back(sum_weight / num_samples);
  }

  return PredictionValues(values, num_types);
}

std::vector<std::pair<double, double>> MultiRegressionPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
