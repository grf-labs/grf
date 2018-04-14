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
                                                             bool use_unweighted_penalty):
        original_data(original_data),
        test_data(test_data),
        lambda(lambda),
        use_unweighted_penalty(use_unweighted_penalty){
};

std::vector<double> LocalLinearPredictionStrategy::predict(size_t sampleID,
                                                           const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                           const Observations& observations) {
  size_t n = observations.get_num_samples();
  size_t num_variables = test_data->get_num_cols();

  Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(n,n);

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it){
    size_t i = it->first;
    double weight = it->second;
    weights(i,i) = weight;
  }

  // generate design matrix X and responses Y as Eigen objects
  Eigen::MatrixXd X(n, num_variables+1);
  Eigen::MatrixXd Y(n, 1);

  // check if we do variable selection. If so, isolate to those variables.

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < num_variables; ++j){
      X(i,j+1) = test_data->get(sampleID, j) - original_data->get(i,j);
    }
    Y(i) = observations.get(Observations::OUTCOME, i);
    X(i, 0) = 1;
  }

  // find ridge regression predictions
  Eigen::MatrixXd M(num_variables+1, num_variables+1);
  M = X.transpose()*weights*X;

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

  Eigen::MatrixXd preds(num_variables+1, 1);
  preds = M.colPivHouseholderQr().solve(X.transpose()*weights*Y);

  return { preds(0) };
}
