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
#include "prediction/LLCausalPredictionStrategy.h"

namespace grf {

LLCausalPredictionStrategy::LLCausalPredictionStrategy(std::vector<double> lambdas,
                                                       bool weight_penalty,
                                                       std::vector<size_t> linear_correction_variables):
        lambdas(lambdas),
        weight_penalty(weight_penalty),
        linear_correction_variables(linear_correction_variables){
};

size_t LLCausalPredictionStrategy::prediction_length() const {
  return lambdas.size();
}

std::vector<double> LLCausalPredictionStrategy::predict(
        size_t sampleID,
        const std::unordered_map<size_t, double>& weights_by_sampleID,
        const Data& train_data,
        const Data& test_data) const {

  // Number of predictor variables to use in local linear regression step
  size_t num_variables = linear_correction_variables.size();

  size_t num_nonzero_weights = weights_by_sampleID.size();
  size_t num_lambdas = lambdas.size();

  // Creating a vector of neighbor weights weights
  // Weights by sample ID contains pairs [sample ID, weight for test point]
  // Weights vec is a vector of weights indexed by their corresponding sample ID
  // For example:
  // weights_by_sampleID = [(0, 0.04), (1, 0.20), (3, 0.01), (4, 0.05), ...]
  // weights_vec = [0.04, 0.20, 0.0, 0.01, 0.05, ....]

  std::vector<size_t> indices(num_nonzero_weights);
  Eigen::MatrixXd weights_vec = Eigen::VectorXd::Zero(num_nonzero_weights);
  {
    size_t i = 0;
    for (const auto& it : weights_by_sampleID) {
      size_t index = it.first;
      double weight = it.second;
      indices[i] = index;
      weights_vec(i) = weight;
      i++;
    }
  }

  // The matrix X consists of differences of linear correction variables from their target.
  // Only observations with nonzero weights need to be filled.
  // For example, if there are K+1 linear correction variables, and m nonzero weights,
  //  then X will be:
  //    1.   (X[0,0] - x[0])   ...   (X[0,K] - x[K]) W[1]
  //    1.   (X[1,0] - x[0])   ...   (X[1,K] - x[K]) W[2]
  //    1.   (X[3,0] - x[0])   ...   (X[3,K] - x[K]) W[3]   # Observation 2 is skipped due to zero weights

  size_t dim_X = 2 * num_variables + 2;
  Eigen::MatrixXd X (num_nonzero_weights, dim_X);
  Eigen::MatrixXd Y (num_nonzero_weights, 1);
  size_t treatment_index = num_variables + 1;

  for (size_t i = 0; i < num_nonzero_weights; ++i) {
    // Index of next neighbor with nonzero weights
    size_t index = indices[i];
    double treatment = train_data.get_treatment(index);

    // Intercept
    X(i, 0) = 1;

    // Regressors
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = linear_correction_variables[j];
      // X - x0 column
      X(i,j+1) = train_data.get(index, current_predictor) -
                    test_data.get(sampleID, current_predictor);
      // (X - x0)*W column
      X(i, treatment_index + j + 1) = X(i, j+1) * treatment;
    }

    // Treatment (just copied)
    X(i, treatment_index) = treatment;

    // Outcome (just copied)
    Y(i) = train_data.get_outcome(index);
  }

  // find ridge regression predictions
  Eigen::MatrixXd M_unpenalized (dim_X, dim_X);
  M_unpenalized.noalias() = X.transpose() * weights_vec.asDiagonal() * X;

  std::vector<double> predictions(num_lambdas);
  Eigen::MatrixXd M;

  for (size_t i = 0; i < num_lambdas; ++i){
    double lambda = lambdas[i];
    M = M_unpenalized;
    if (!weight_penalty) {
      double normalization = M_unpenalized.trace() / dim_X;

      // standard ridge penalty
      for (size_t j = 1; j < dim_X; ++j){
        if (j != treatment_index){
          M(j, j) += lambda * normalization;
        }
      }
    } else {
      // covariance ridge penalty
      for (size_t j = 1; j < dim_X; ++j){
        if (j != treatment_index){
          M(j, j) += lambda * M(j, j);
        }
      }
    }

    Eigen::MatrixXd local_coefficients = M.ldlt().solve(X.transpose()*weights_vec.asDiagonal()*Y);

    // We're only interested in the coefficient associated with the treatment variable
    predictions[i] = local_coefficients(treatment_index);
  }

  return predictions;
}

std::vector<double> LLCausalPredictionStrategy::compute_variance(
        size_t sampleID,
        const std::vector<std::vector<size_t>>& samples_by_tree,
        const std::unordered_map<size_t, double>& weights_by_sampleID,
        const Data& train_data,
        const Data& test_data,
        size_t ci_group_size) const {

  double lambda = lambdas[0];

  size_t num_variables = linear_correction_variables.size();
  size_t num_nonzero_weights = weights_by_sampleID.size();

  std::vector<size_t> sample_index_map(train_data.get_num_rows());
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

  size_t dim_X = 2 * num_variables + 2;
  Eigen::MatrixXd X (num_nonzero_weights, dim_X);
  Eigen::MatrixXd Y (num_nonzero_weights, 1);
  size_t treatment_index = num_variables + 1;

  for (size_t i = 0; i < num_nonzero_weights; ++i) {
    // Index of next neighbor with nonzero weights

    size_t index = indices[i];
    double treatment = train_data.get_treatment(index);

    // Intercept
    X(i, 0) = 1;

    // Regressors
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = linear_correction_variables[j];
      // X - x0 column
      X(i,j+1) = train_data.get(index, current_predictor) -
                 test_data.get(sampleID, current_predictor);

      // (X - x0)*W column
      X(i, treatment_index + j + 1) = X(i, j+1) * treatment;
    }

    // Treatment (just copied)
    X(i, treatment_index) = treatment;

    // Outcome (just copied)
    Y(i) = train_data.get_outcome(index);
  }

  // find ridge regression predictions
  Eigen::MatrixXd M_unpenalized (dim_X, dim_X);
  M_unpenalized.noalias() = X.transpose() * weights_vec.asDiagonal() * X;

  Eigen::MatrixXd M;
  M = M_unpenalized;
  if (!weight_penalty) {
    double normalization = M_unpenalized.trace() / dim_X;

    // standard ridge penalty
    for(size_t j = 1; j < dim_X; ++j){
      if(j != treatment_index){
        M(j, j) += lambda * normalization;
      }
    }
  } else {
    // covariance ridge penalty
    for(size_t j = 1; j < dim_X; ++j){
      if(j != treatment_index){
        M(j, j) += lambda * M(j, j);
      }
    }
  }

  Eigen::VectorXd theta = M.ldlt().solve(X.transpose()*weights_vec.asDiagonal()*Y);

  Eigen::VectorXd e_trt = Eigen::VectorXd::Zero(dim_X);
  e_trt(treatment_index) = 1.0;
  Eigen::VectorXd zeta = M.ldlt().solve(e_trt);

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
      for (size_t sample : samples_by_tree[b]) {
        psi_1 += pseudo_residual(sample_index_map[sample]);
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

} // namespace grf
