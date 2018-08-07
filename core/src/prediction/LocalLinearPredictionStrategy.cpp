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
      X(i,j+1) = test_data->get(sampleID, current_predictor)
              - original_data->get(indices[i],current_predictor);
    }
    Y(i) = observations.get(Observations::OUTCOME, indices[i]);
    X(i, 0) = 1;
  }

  // find ridge regression predictions
  Eigen::MatrixXd M (num_variables+1, num_variables+1);
  M.noalias() = X.transpose()*weights_vec.asDiagonal()*X;

  size_t num_lambdas = lambdas.size();
  std::vector<double> predictions(num_lambdas);

  for( size_t i = 0; i < num_lambdas; ++i){
    double lambda = lambdas[i];
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
