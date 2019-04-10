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
#include <string>
#include <vector>
#include <iostream>
#include "Eigen/Dense"
#include "commons/utility.h"
#include "commons/Data.h"
#include "prediction/LocalLinearPredictionStrategy.h"

size_t LocalLinearPredictionStrategy::prediction_length() {
  return lambdas.size();
}

LocalLinearPredictionStrategy::LocalLinearPredictionStrategy(std::vector<double> lambdas,
                                                             bool weight_penalty,
                                                             std::vector<size_t> linear_correction_variables):
        lambdas(lambdas),
        weight_penalty(weight_penalty),
        linear_correction_variables(linear_correction_variables){
};

std::vector<double> LocalLinearPredictionStrategy::predict(
    size_t sampleID,
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    const Data* train_data,
    const Data* data) {
  size_t num_variables = linear_correction_variables.size();
  size_t num_nonzero_weights = weights_by_sampleID.size();

  std::vector<size_t> indices(num_nonzero_weights);
  Eigen::MatrixXd weights_vec = Eigen::VectorXd::Zero(num_nonzero_weights);
  {
    size_t i = 0;
    for (auto& it : weights_by_sampleID) {
      size_t index = it.first;
      double weight = it.second;
      indices[i] = index;
      weights_vec(i) = weight;
      i++;
    }
  }

  Eigen::MatrixXd X (num_nonzero_weights, num_variables+1);
  Eigen::MatrixXd Y (num_nonzero_weights, 1);
  for (size_t i = 0; i < num_nonzero_weights; ++i) {
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = linear_correction_variables[j];
      X(i,j+1) = train_data->get(indices[i],current_predictor)
                 - data->get(sampleID, current_predictor);
    }
    Y(i) = train_data->get_outcome(indices[i]);
    X(i, 0) = 1;
  }

  // find ridge regression predictions
  Eigen::MatrixXd M_unpenalized(num_variables+1, num_variables+1);
  M_unpenalized.noalias() = X.transpose()*weights_vec.asDiagonal()*X;

  Eigen::MatrixXd M;
  size_t num_lambdas = lambdas.size();
  std::vector<double> predictions(num_lambdas);

  for (size_t i = 0; i < num_lambdas; ++i){
    double lambda = lambdas[i];
    M = M_unpenalized;

    if (!weight_penalty) {
      // standard ridge penalty
      double normalization = M.trace() / (num_variables + 1);
      for (size_t j = 1; j < num_variables + 1; ++j){
        M(j,j) += lambda * normalization;
      }
    } else {
      // covariance ridge penalty
      for (size_t j = 1; j < num_variables+1; ++j){
        M(j,j) += lambda * M(j,j); // note that the weights are already normalized
      }
    }

    Eigen::MatrixXd local_coefficients = M.ldlt().solve(X.transpose()*weights_vec.asDiagonal()*Y);
    predictions[i] =  local_coefficients(0);
  }

  return predictions;
}

std::vector<double> LocalLinearPredictionStrategy::compute_variance(
    size_t sampleID,
    std::vector<std::vector<size_t>> samples_by_tree,
    std::unordered_map<size_t, double> weights_by_sampleID,
    const Data* train_data,
    const Data* data,
    size_t ci_group_size) {

  double lambda = lambdas[0];

  size_t num_variables = linear_correction_variables.size();
  size_t num_nonzero_weights = weights_by_sampleID.size();

  std::vector<size_t> sample_index_map(train_data->get_num_rows());
  std::vector<size_t> indices(num_nonzero_weights);

  Eigen::MatrixXd weights_vec = Eigen::VectorXd::Zero(num_nonzero_weights);
    {
      size_t i = 0;
      for (auto& it : weights_by_sampleID) {
        size_t index = it.first;
        double weight = it.second;
        indices[i] = index;
        sample_index_map[index] = i;
        weights_vec(i) = weight;
        i++;
      }
  }

  Eigen::MatrixXd X (num_nonzero_weights, num_variables+1);
  Eigen::MatrixXd Y (num_nonzero_weights, 1);
  for (size_t i = 0; i < num_nonzero_weights; ++i) {
    X(i, 0) = 1;
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = linear_correction_variables[j];
      X(i,j+1) = train_data->get(indices[i],current_predictor)
                 - data->get(sampleID, current_predictor);
    }
    Y(i) = train_data->get_outcome(indices[i]);
  }

  // find ridge regression predictions
  Eigen::MatrixXd M (num_variables+1, num_variables+1);
  M.noalias() = X.transpose()*weights_vec.asDiagonal()*X;

  if (!weight_penalty) {
    double normalization = M.trace() / (num_variables + 1);
    for (size_t i = 1; i < num_variables + 1; ++i){
      M(i,i) += lambda * normalization;
    }
  } else {
    for (size_t i = 1; i < num_variables+1; ++i){
      M(i,i) += lambda * M(i,i);
    }
  }

  Eigen::VectorXd theta = M.ldlt().solve(X.transpose()*weights_vec.asDiagonal()*Y);

  Eigen::VectorXd e_one = Eigen::VectorXd::Zero(num_variables+1);
  e_one(0) = 1.0;
  Eigen::VectorXd zeta = M.ldlt().solve(e_one);

  Eigen::VectorXd X_times_zeta = X * zeta;
  Eigen::VectorXd local_prediction = X * theta;
  Eigen::VectorXd pseudo_residual = Eigen::VectorXd::Zero(num_nonzero_weights);

  for (size_t i = 0; i < num_nonzero_weights; i++) {
    pseudo_residual(i) = X_times_zeta(i) * (Y(i) - local_prediction(i));
  }

  double num_good_groups = 0;
  double psi_squared = 0;
  double psi_grouped_squared = 0;

  double avg_score = 0;

  for (size_t group = 0; group < samples_by_tree.size() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (samples_by_tree[group * ci_group_size + j].size() == 0) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_psi = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {
      size_t b = group * ci_group_size + j;
      double psi_1 = 0;
      for(size_t k = 0; k < samples_by_tree[b].size(); ++ k){
        psi_1 += pseudo_residual(sample_index_map[samples_by_tree[b][k]]);
      }
      psi_1 /= samples_by_tree[b].size();
      psi_squared += psi_1 * psi_1;
      group_psi += psi_1;
    }

    group_psi /= ci_group_size;
    psi_grouped_squared += group_psi * group_psi;

    avg_score += group_psi;
  }

  avg_score /= num_good_groups;

  double var_between = psi_grouped_squared / num_good_groups - avg_score * avg_score;
  double var_total = psi_squared / (num_good_groups * ci_group_size) - avg_score * avg_score;

  // This is the amount by which var_between is inflated due to using small groups
  double group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an objective
  // Bayes analysis of variance instead to avoid negative values.
  double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);

  return { var_debiased };
}
