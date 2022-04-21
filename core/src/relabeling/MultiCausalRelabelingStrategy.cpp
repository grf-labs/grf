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

#include "commons/utility.h"
#include "relabeling/MultiCausalRelabelingStrategy.h"

namespace grf {

MultiCausalRelabelingStrategy::MultiCausalRelabelingStrategy(size_t response_length,
                                                             const std::vector<double>& gradient_weights) {
  this->response_length = response_length;
  if (gradient_weights.empty()) {
    this->gradient_weights = std::vector<double> (response_length, 1.0);
  } else if (gradient_weights.size() != response_length) {
    throw std::runtime_error("Optional gradient weights vector must be same length as response_length.");
  } else {
    this->gradient_weights = gradient_weights;
  }
}

bool MultiCausalRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    Eigen::ArrayXXd& responses_by_sample) const {

  // Prepare the relevant averages.
  size_t num_samples = samples.size();
  size_t num_treatments = data.get_num_treatments();
  size_t num_outcomes = data.get_num_outcomes();
  if (num_samples <= num_treatments) {
    return true;
  }

  Eigen::MatrixXd Y_centered = Eigen::MatrixXd(num_samples, num_outcomes);
  Eigen::MatrixXd W_centered = Eigen::MatrixXd(num_samples, num_treatments);
  Eigen::VectorXd weights = Eigen::VectorXd(num_samples);
  Eigen::VectorXd Y_mean = Eigen::VectorXd::Zero(num_outcomes);
  Eigen::VectorXd W_mean = Eigen::VectorXd::Zero(num_treatments);
  double sum_weight = 0;
  for (size_t i = 0; i < num_samples; i++) {
    size_t sample = samples[i];
    double weight = data.get_weight(sample);
    Eigen::VectorXd outcome = data.get_outcomes(sample);
    Eigen::VectorXd treatment = data.get_treatments(sample);
    Y_centered.row(i) = outcome;
    W_centered.row(i) = treatment;
    weights(i) = weight;
    Y_mean += weight * outcome;
    W_mean += weight * treatment;
    sum_weight += weight;
  }
  Y_mean /= sum_weight;
  W_mean /= sum_weight;
  Y_centered.rowwise() -= Y_mean.transpose();
  W_centered.rowwise() -= W_mean.transpose();

  if (std::abs(sum_weight) <= 1e-16) {
    return true;
  }

  Eigen::MatrixXd WW_bar = W_centered.transpose() * weights.asDiagonal() * W_centered; // [num_treatments X num_treatments]
  // Calculate the treatment effect.
  // This condition number check works fine in practice - there may be more robust ways.
  if (equal_doubles(WW_bar.determinant(), 0.0, 1.0e-10)) {
    return true;
  }

  Eigen::MatrixXd A_p_inv = WW_bar.inverse();
  Eigen::MatrixXd beta = A_p_inv * W_centered.transpose() * weights.asDiagonal() * Y_centered; // [num_treatments X num_outcomes]

  Eigen::MatrixXd rho_weight = W_centered * A_p_inv.transpose(); // [num_samples X num_treatments]
  Eigen::MatrixXd residual = Y_centered - W_centered * beta; // [num_samples X num_outcomes]

  // Create the new outcomes, eq (20) in https://arxiv.org/pdf/1610.01271.pdf
  // `responses_by_sample(sample_i, )` is a `num_treatments*num_outcomes`-sized vector.
  for (size_t i = 0; i < num_samples; i++) {
    size_t sample = samples[i];
    size_t j = 0;
    for (size_t outcome = 0; outcome < num_outcomes; outcome++) {
      for (size_t treatment = 0; treatment < num_treatments; treatment++) {
        responses_by_sample(sample, j) = rho_weight(i, treatment) * residual(i, outcome) * gradient_weights[j];
        j++;
      }
    }
  }
  return false;
}

size_t MultiCausalRelabelingStrategy::get_response_length() const {
  return response_length;
}

} // namespace grf
