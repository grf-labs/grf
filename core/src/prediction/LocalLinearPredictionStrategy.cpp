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
#include "commons/Observations.h"
#include "prediction/LocalLinearPredictionStrategy.h"


const size_t LocalLinearPredictionStrategy::OUTCOME = 0;

size_t LocalLinearPredictionStrategy::prediction_length() {
  return lambdas.size();
}

LocalLinearPredictionStrategy::LocalLinearPredictionStrategy(const Data* original_data,
                                                             const Data* test_data,
                                                             std::vector<double> lambdas,
                                                             bool use_unweighted_penalty,
                                                             std::vector<size_t> linear_correction_variables):
        original_data(original_data),
        test_data(test_data),
        lambdas(lambdas),
        use_unweighted_penalty(use_unweighted_penalty),
        linear_correction_variables(linear_correction_variables){
};

std::vector<double> LocalLinearPredictionStrategy::predict(
    size_t sampleID,
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    const Observations& observations) {
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
      X(i,j+1) = original_data->get(indices[i],current_predictor)
                 - test_data->get(sampleID, current_predictor);
    }
    Y(i) = observations.get(Observations::OUTCOME, indices[i]);
    X(i, 0) = 1;
  }

  // find ridge regression predictions
  Eigen::MatrixXd M (num_variables+1, num_variables+1);
  M.noalias() = X.transpose()*weights_vec.asDiagonal()*X;
  Eigen::MatrixXd M_stored = M;

  size_t num_lambdas = lambdas.size();
  std::vector<double> predictions(num_lambdas);

  for( size_t i = 0; i < num_lambdas; ++i){
    double lambda = lambdas[i];

    M = M_stored;
    if (use_unweighted_penalty) {
      // standard ridge penalty
      double normalization = M.trace() / (num_variables + 1);
      for (size_t i = 1; i < num_variables + 1; ++i){
        M(i,i) += lambda * normalization;
      }
    } else {
      // covariance ridge penalty
      for (size_t i = 1; i < num_variables+1; ++i){
        M(i,i) += lambda * M(i,i); // note that the weights are already normalized
      }
    }

    Eigen::MatrixXd preds = M.ldlt().solve(X.transpose()*weights_vec.asDiagonal()*Y);
    predictions[i] =  preds(0);
  }

  return predictions;
}

Eigen::MatrixXd LocalLinearPredictionStrategy::find_M(
        std::unordered_map<size_t, double> weights_by_sampleID,
        size_t sampleID,
        double lambda){

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
  for (size_t i = 0; i < num_nonzero_weights; ++i) {
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = linear_correction_variables[j];
      X(i,j+1) = original_data->get(indices[i],current_predictor)
                 - test_data->get(sampleID, current_predictor);
    }
    X(i, 0) = 1;
  }

  // find ridge regression predictions
  Eigen::MatrixXd M (num_variables+1, num_variables+1);
  M.noalias() = X.transpose()*weights_vec.asDiagonal()*X;

  if (use_unweighted_penalty) {
    // standard ridge penalty
    double normalization = M.trace() / (num_variables + 1);
    for (size_t i = 1; i < num_variables + 1; ++i){
      M(i,i) += lambda * normalization;
    }
  } else {
    // covariance ridge penalty
    for (size_t i = 1; i < num_variables+1; ++i){
      M(i,i) += lambda * M(i,i);
    }
  }
  return M;
}

std::vector<double> LocalLinearPredictionStrategy::compute_variance(
        const PredictionValues& leaf_values,
        uint ci_group_size,
        size_t sampleID,
        std::unordered_map<size_t, double> weights_by_sampleID,
        double prediction_sample){

  size_t num_variables = linear_correction_variables.size();

  double lambda = lambdas[0];
  Eigen::MatrixXd M_total = find_M(weights_by_sampleID, sampleID, lambda);
  Eigen::VectorXd e_one = Eigen::VectorXd::Zero(num_variables+1);
  e_one(0) = 1.0;
  Eigen::MatrixXd M_inverse = M_total.ldlt().solve(e_one);

  Eigen::VectorXd X_sample (1, num_variables);
  for(size_t index = 0; index < num_variables; ++ index){
    size_t current_predictor = linear_correction_variables[index];
    X_sample(index) = test_data->get(sampleID, current_predictor);
  }

  Eigen::MatrixXd M_X = M_inverse * X_sample;
  double factor = M_X(0);

  double num_good_groups = 0;
  double psi_squared = 0;
  double psi_grouped_squared = 0;

  for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_values.empty(group * ci_group_size + j)) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_psi = 0;
    double psi_squared = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {
      size_t i = group * ci_group_size + j;

      double prediction_group = leaf_values.get(i, OUTCOME);
      double psi_1 = factor * (prediction_sample - prediction_group);

      psi_squared += psi_1 * psi_1;
      group_psi += psi_1;
    }

    group_psi /= ci_group_size;
    psi_grouped_squared += group_psi * group_psi;
  }

  double var_between = psi_grouped_squared / num_good_groups;
  double var_total = psi_squared / (num_good_groups * ci_group_size);

  // This is the amount by which var_between is inflated due to using small groups
  double group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an objective
  // Bayes analysis of variance instead to avoid negative values.
  double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);

  return { var_debiased };
}

std::vector<double> LocalLinearPredictionStrategy::compute_debiased_error(
        size_t sample,
        const PredictionValues& leaf_values,
        const Observations& observations) {
  double outcome = observations.get(Observations::OUTCOME, sample);

  double error = outcome;
  double mse = error * error;

  double bias = 0.0;
  size_t num_trees = 0;
  for (size_t n = 0; n < leaf_values.get_num_nodes(); n++) {
    if (leaf_values.empty(n)) {
      continue;
    }

    double tree_variance = leaf_values.get(n, OUTCOME);
    bias += tree_variance * tree_variance;
    num_trees++;
  }

  if (num_trees <= 1) {
    return { NAN };
  }

  bias /= num_trees * (num_trees - 1);
  return { mse - bias };
}

size_t LocalLinearPredictionStrategy::prediction_value_length() {
  return 1;
}
