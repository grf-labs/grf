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

#include "relabeling/LLRegressionRelabelingStrategy.h"

namespace grf {

LLRegressionRelabelingStrategy::LLRegressionRelabelingStrategy(double split_lambda,
                                                               bool weight_penalty,
                                                               const std::vector<double>& overall_beta,
                                                               size_t ll_split_cutoff,
                                                               std::vector<size_t> ll_split_variables):
  split_lambda(split_lambda),
  weight_penalty(weight_penalty),
  overall_beta(overall_beta),
  ll_split_cutoff(ll_split_cutoff),
  ll_split_variables(ll_split_variables){
};

bool LLRegressionRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    Eigen::ArrayXXd& responses_by_sample) const {

  size_t num_variables = ll_split_variables.size();
  size_t num_data_points = samples.size();

  Eigen::MatrixXd X (num_data_points, num_variables+1);
  Eigen::MatrixXd Y (num_data_points, 1);
  for (size_t i = 0; i < num_data_points; ++i) {
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = ll_split_variables[j];
      X(i, j + 1) = data.get(samples[i],current_predictor);
    }
    Y(i) = data.get_outcome(samples[i]);
    X(i, 0) = 1;
  }

  Eigen::MatrixXd leaf_predictions (num_data_points, 1);

  if (num_data_points < ll_split_cutoff) {
    // use overall beta for ridge predictions

    Eigen::MatrixXd eigen_beta (num_variables + 1, 1);
    for(size_t j = 0; j < num_variables + 1; ++j){
      eigen_beta(j) = overall_beta[j];
    }
    leaf_predictions = X * eigen_beta;
  } else {
    // find ridge regression predictions
    Eigen::MatrixXd M(num_variables + 1, num_variables + 1);
    M.noalias() = X.transpose() * X;

    if (!weight_penalty) {
      // standard ridge penalty
      double normalization = M.trace() / (num_variables + 1);
      for (size_t j = 1; j < num_variables + 1; ++j){
        M(j, j) += split_lambda * normalization;
      }
    } else {
      // covariance ridge penalty
      for (size_t j = 1; j < num_variables + 1; ++j){
        M(j, j) += split_lambda * M(j,j); // note that the weights are already normalized
      }
    }

    Eigen::MatrixXd local_coefficients = M.ldlt().solve(X.transpose() * Y);
    leaf_predictions = X * local_coefficients;
  }

  size_t i = 0;
  for (size_t sample : samples) {
      double prediction_sample = leaf_predictions(i);
      double residual = prediction_sample - data.get_outcome(sample);
      responses_by_sample(sample, 0) = residual;
      i++;
  }
    return false;
  }

} // namespace grf
