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
Rcpp::List causal_train(Rcpp::NumericMatrix train_matrix,
                        Eigen::SparseMatrix<double> sparse_train_matrix,
                        size_t outcome_index,
                        size_t treatment_index,
                        size_t sample_weight_index,
                        bool use_sample_weights,
                        unsigned int mtry,
                        unsigned int num_trees,
                        unsigned int min_node_size,
                        double sample_fraction,
                        bool honesty,
                        double honesty_fraction,
                        size_t ci_group_size,
                        double reduced_form_weight,
                        double alpha,
                        double imbalance_penalty,
                        bool stabilize_splits,
                        std::vector<size_t> clusters,
                        unsigned int samples_per_cluster,
                        bool compute_oob_predictions,
                        unsigned int num_threads,
                        unsigned int seed) {
  ForestTrainer trainer = ForestTrainers::instrumental_trainer(reduced_form_weight, stabilize_splits);

  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_treatment_index(treatment_index - 1);
  data->set_instrument_index(treatment_index - 1);
  if(use_sample_weights) {
      data->set_weight_index(sample_weight_index - 1);
  }
  data->sort();

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
                        honesty_fraction, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads);
    predictions = predictor.predict_oob(forest, data, false);
  }

  delete data;
  return RcppUtilities::create_forest_object(forest, predictions);
}


// [[Rcpp::export]]
Rcpp::List causal_predict(Rcpp::List forest_object,
                          Rcpp::NumericMatrix train_matrix,
                          Eigen::SparseMatrix<double> sparse_train_matrix,
                          size_t outcome_index,
                          size_t treatment_index,
                          Rcpp::NumericMatrix test_matrix,
                          Eigen::SparseMatrix<double> sparse_test_matrix,
                          unsigned int num_threads,
                          bool estimate_variance) {
  Data* train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  train_data->set_outcome_index(outcome_index - 1);
  train_data->set_treatment_index(treatment_index - 1);
  train_data->set_instrument_index(treatment_index - 1);
  Data* data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

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
                              unsigned int num_threads,
                              bool estimate_variance) {
  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_treatment_index(treatment_index - 1);
  data->set_instrument_index(treatment_index - 1);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

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
                                 size_t outcome_index,
                                 size_t treatment_index,
                                 std::vector<double> lambdas,
                                 bool use_weighted_penalty,
                                 std::vector<size_t> linear_correction_variables,
                                 unsigned int num_threads,
                                 bool estimate_variance) {
  Data* test_data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Data* train_data = RcppUtilities::convert_data(training_data, sparse_training_data);

  train_data->set_outcome_index(outcome_index - 1);
  train_data->set_treatment_index(treatment_index - 1);
  train_data->set_instrument_index(treatment_index - 1);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest);

  ForestPredictor predictor = ForestPredictors::ll_causal_predictor(num_threads, lambdas, use_weighted_penalty,
                                                                 linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, train_data, test_data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete train_data;
  delete test_data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List ll_causal_predict_oob(Rcpp::List forest,
                                     Rcpp::NumericMatrix input_data,
                                     Eigen::SparseMatrix<double> sparse_input_data,
                                     size_t outcome_index,
                                     size_t treatment_index,
                                     std::vector<double> lambdas,
                                     bool use_weighted_penalty,
                                     std::vector<size_t> linear_correction_variables,
                                     unsigned int num_threads,
                                     bool estimate_variance) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);

  data->set_outcome_index(outcome_index - 1);
  data->set_treatment_index(treatment_index - 1);
  data->set_instrument_index(treatment_index - 1);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest);

  ForestPredictor predictor = ForestPredictors::ll_causal_predictor(num_threads, lambdas, use_weighted_penalty,
                                                                 linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}
