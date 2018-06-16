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
  return 1;
}

LocalLinearPredictionStrategy::LocalLinearPredictionStrategy(const Data* original_data,
                                                             const Data* test_data,
                                                             double lambda,
                                                             bool use_unweighted_penalty,
                                                             std::vector<size_t> linear_correction_variables):
        original_data(original_data),
        test_data(test_data),
        lambda(lambda),
        use_unweighted_penalty(use_unweighted_penalty),
        linear_correction_variables(linear_correction_variables){
};

std::vector<double> LocalLinearPredictionStrategy::predict(size_t sampleID,
                                                           const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                           const Observations& observations) {
  size_t num_variables = linear_correction_variables.size();
  size_t nNZ = 0;
  for (auto const & it : weights_by_sampleID) {
    double weight = it.second;
    if (weight > 0) {
      nNZ++;
    }
  }
  std::vector<size_t> indices;
  Eigen::MatrixXd weightsVec = Eigen::VectorXd::Zero(nNZ);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(nNZ, num_variables+1);
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(nNZ, 1);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(num_variables+1, num_variables+1);
  Eigen::MatrixXd preds = Eigen::MatrixXd::Zero(num_variables+1, 1);

  size_t curIdx = 0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it){
    size_t i = it->first;
    double weight = it->second;
    if (weight > 0) {
      indices.push_back(i);
      weightsVec(curIdx) = weight;
      curIdx++;
    }
  }

  for (size_t i = 0; i < nNZ; ++i) {
    for (size_t j = 0; j < num_variables; ++j){
      size_t current_predictor = linear_correction_variables[j];
      X(i,j+1) = test_data->get(sampleID, current_predictor)
              - original_data->get(indices[i],current_predictor);
    }
    Y(i) = observations.get(Observations::OUTCOME, indices[i]);
    X(i, 0) = 1;
  }

  // find ridge regression predictions
  M.noalias() = X.transpose()*weightsVec.asDiagonal()*X;

  if (use_unweighted_penalty) {
    // standard ridge penalty
    double additional_regularization = M.trace() / (num_variables + 1);
    for (size_t i = 1; i < num_variables + 1; ++i){
      M(i,i) += lambda * additional_regularization;
    }
  } else {
    // covariance ridge penalty
    for (size_t i = 1; i < num_variables+1; ++i){
      M(i,i) += lambda * M(i,i); // note that the weights are already normalized
    }
  }

  preds = M.ldlt().solve(X.transpose()*weightsVec.asDiagonal()*Y);

  return { preds(0) };
}
