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

#include <cmath>
#include <vector>

#include "commons/Data.h"
#include "commons/utility.h"
#include "prediction/MultiCausalPredictionStrategy.h"

namespace grf {

MultiCausalPredictionStrategy::MultiCausalPredictionStrategy(size_t num_treatments):
  num_treatments(num_treatments) {}

size_t MultiCausalPredictionStrategy::prediction_length() const {
    return num_treatments;
}

std::vector<double> MultiCausalPredictionStrategy::predict(const std::vector<double>& average) const {
  return average;
}

std::vector<double> MultiCausalPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {
  return { 0.0 };
}

size_t MultiCausalPredictionStrategy::prediction_value_length() const {
  return 1;
}

PredictionValues MultiCausalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  // Compute the sample weighted least squares solution to Y = W * beta + epsilon
  for (size_t i = 0; i < leaf_samples.size(); ++i) {
    const std::vector<size_t>& samples = leaf_samples[i];
    size_t num_samples = samples.size();
    if (num_samples <= num_treatments) {
      continue;
    }

    Eigen::VectorXd Y_centered = Eigen::VectorXd(num_samples);
    Eigen::MatrixXd W_centered = Eigen::MatrixXd(num_samples, num_treatments);
    Eigen::VectorXd weights = Eigen::VectorXd(num_samples);
    double Y_mean = 0;
    Eigen::VectorXd W_mean = Eigen::VectorXd::Zero(num_treatments);
    double sum_weight = 0;
    for (size_t j = 0; j < num_samples; j++) {
      size_t sample = samples[j];
      double weight = data.get_weight(sample);
      double outcome = data.get_outcome(sample);
      Eigen::VectorXd treatment = data.get_treatments(sample);
      Y_centered(j) = outcome;
      W_centered.row(j) = treatment;
      weights(j) = weight;
      Y_mean += weight * outcome;
      W_mean += weight * treatment;
      sum_weight += weight;
    }
    Y_mean = Y_mean / sum_weight;
    W_mean = W_mean / sum_weight;
    Y_centered.array() -= Y_mean;
    W_centered.rowwise() -= W_mean.transpose();

    // if total weight is very small, treat the leaf as empty
    if (std::abs(sum_weight) <= 1e-16) {
      continue;
    }

    // The following "condition" check appear to work reasonably well in practice
    if (equal_doubles((W_centered.transpose() * weights.asDiagonal() * W_centered).determinant(), 0.0, 1.0e-2)) {
      continue;
    }

    // Solving W beta = Y is the same as W^T W beta = W^T Y and provides a moderate speedup over
    // the matrix inverse. We are in a leaf node with a limited number of samples so this is not a very
    // performance sensitive concern.
    weights = weights.array().sqrt().matrix();
    Eigen::VectorXd beta = ((weights.asDiagonal() * W_centered).transpose() * weights.asDiagonal() * W_centered).
      ldlt().
      solve((weights.asDiagonal() * W_centered).transpose() * weights.asDiagonal() * Y_centered);

    values[i] = std::vector<double> (beta.data(), beta.data() + num_treatments);
  }

  return PredictionValues(values, num_treatments);
}

std::vector<std::pair<double, double>> MultiCausalPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
