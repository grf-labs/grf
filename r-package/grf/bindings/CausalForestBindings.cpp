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

#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List causal_predict(Rcpp::List forest_object,
                          Rcpp::NumericMatrix train_matrix,
                          Eigen::SparseMatrix<double> sparse_train_matrix,
                          size_t outcome_index,
                          size_t treatment_index,
                          size_t instrument_index,
                          Rcpp::NumericMatrix test_matrix,
                          Eigen::SparseMatrix<double> sparse_test_matrix,
                          unsigned int num_threads,
                          bool estimate_variance) {
  Data* train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  train_data->set_outcome_index(outcome_index - 1);
  train_data->set_treatment_index(treatment_index - 1);
  train_data->set_instrument_index(instrument_index - 1);
  Data* data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);

  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, train_data, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete train_data;
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List causal_predict_oob(Rcpp::List forest_object,
                              Rcpp::NumericMatrix train_matrix,
                              Eigen::SparseMatrix<double> sparse_train_matrix,
                              size_t outcome_index,
                              size_t treatment_index,
                              size_t instrument_index,
                              unsigned int num_threads,
                              bool estimate_variance) {
  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_treatment_index(treatment_index - 1);
  data->set_instrument_index(instrument_index - 1);

  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List ll_causal_predict(Rcpp::List forest,
                                 Rcpp::NumericMatrix input_data,
                                 Rcpp::NumericMatrix training_data,
                                 Eigen::SparseMatrix<double> sparse_input_data,
                                 Eigen::SparseMatrix<double> sparse_training_data,
                                 std::vector<double> lambdas,
                                 bool use_unweighted_penalty,
                                 std::vector<size_t> linear_correction_variables,
                                 unsigned int num_threads) {
  Data *test_data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Data *train_data = RcppUtilities::convert_data(training_data, sparse_training_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::ll_causal_predictor(num_threads, lambdas, use_unweighted_penalty,
                                                                 linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, train_data, test_data, false);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete train_data;
  delete test_data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List ll_causal_predict_oob(Rcpp::List forest,
                                     Rcpp::NumericMatrix input_data,
                                     Eigen::SparseMatrix<double> sparse_input_data,
                                     std::vector<double> lambdas,
                                     bool use_unweighted_penalty,
                                     std::vector<size_t> linear_correction_variables,
                                     unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::ll_causal_predictor(num_threads, lambdas, use_unweighted_penalty,
                                                                 linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data, false);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}
