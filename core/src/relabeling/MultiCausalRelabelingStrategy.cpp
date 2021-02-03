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

bool MultiCausalRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    Eigen::ArrayXXd& responses_by_sample) const {

  // Prepare the relevant averages.
  size_t num_samples = samples.size();
  size_t num_treatments = data.get_num_treatments();
  if (num_samples <= num_treatments) {
    return true;
  }

  Eigen::VectorXd Y_centered = Eigen::VectorXd(num_samples);
  Eigen::MatrixXd W_centered = Eigen::MatrixXd(num_samples, num_treatments);
  Eigen::VectorXd weights = Eigen::VectorXd(num_samples);
  double Y_mean = 0;
  Eigen::VectorXd W_mean = Eigen::VectorXd::Zero(num_treatments);
  double sum_weight = 0;
  for (size_t i = 0; i < num_samples; i++) {
    size_t sample = samples[i];
    double weight = data.get_weight(sample);
    double outcome = data.get_outcome(sample);
    Eigen::VectorXd treatment = data.get_treatments(sample);
    Y_centered(i) = outcome;
    W_centered.row(i) = treatment;
    weights(i) = weight;
    Y_mean += weight * outcome;
    W_mean += weight * treatment;
    sum_weight += weight;
  }
  Y_mean = Y_mean / sum_weight;
  W_mean = W_mean / sum_weight;
  Y_centered.array() -= Y_mean;
  W_centered.rowwise() -= W_mean.transpose();

  if (std::abs(sum_weight) <= 1e-16) {
    return true;
  }

  // Calculate the treatment effect.
  // This condition number check works fine in practice - there may be more robust ways.
  if (equal_doubles((W_centered.transpose() * weights.asDiagonal() * W_centered).determinant(), 0.0, 1.0e-10)) {
    return true;
  }

  Eigen::MatrixXd A_p_inv = (W_centered.transpose() * weights.asDiagonal() * W_centered).inverse();
  Eigen::VectorXd beta = A_p_inv * W_centered.transpose() * weights.asDiagonal() * Y_centered;

  Eigen::ArrayXXd rho = (W_centered * A_p_inv.transpose()).array().colwise() * (Y_centered - W_centered * beta).array();

  // Create the new outcomes.
  for (size_t i = 0; i < num_samples; i++) {
    size_t sample = samples[i];
    responses_by_sample.row(sample) = rho.row(i);
  }
  return false;
}

} // namespace grf
