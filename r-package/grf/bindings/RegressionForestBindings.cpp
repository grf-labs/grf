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

#include <Rcpp.h>
#include <vector>

#include "commons/globals.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List regression_train(const Rcpp::NumericMatrix& train_matrix,
                            size_t outcome_index,
                            size_t sample_weight_index,
                            bool use_sample_weights,
                            unsigned int mtry,
                            unsigned int num_trees,
                            unsigned int min_node_size,
                            double sample_fraction,
                            bool honesty,
                            double honesty_fraction,
                            bool honesty_prune_leaves,
                            size_t ci_group_size,
                            double alpha,
                            double imbalance_penalty,
                            std::vector<size_t> clusters,
                            unsigned int samples_per_cluster,
                            bool compute_oob_predictions,
                            unsigned int num_threads,
                            unsigned int seed) {
  ForestTrainer trainer = regression_trainer();

  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_outcome_index(outcome_index);
  if(use_sample_weights) {
      data.set_weight_index(sample_weight_index);
  }

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = regression_predictor(num_threads);
    predictions = predictor.predict_oob(forest, data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List regression_predict(Rcpp::List forest_object,
                              const Rcpp::NumericMatrix& train_matrix,
                              size_t outcome_index,
                              const Rcpp::NumericMatrix& test_matrix,
                              unsigned int num_threads,
                              unsigned int estimate_variance) {
  Data train_data = RcppUtilities::convert_data(train_matrix);
  train_data.set_outcome_index(outcome_index);

  Data data = RcppUtilities::convert_data(test_matrix);
  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = regression_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, train_data, data, estimate_variance);

  return RcppUtilities::create_prediction_object(predictions);
}

// [[Rcpp::export]]
Rcpp::List regression_predict_oob(Rcpp::List forest_object,
                                  const Rcpp::NumericMatrix& train_matrix,
                                  size_t outcome_index,
                                  unsigned int num_threads,
                                  bool estimate_variance) {
  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_outcome_index(outcome_index);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = regression_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, estimate_variance);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  return result;
}

// [[Rcpp::export]]
Rcpp::List ll_regression_train(const Rcpp::NumericMatrix& train_matrix,
                            size_t outcome_index,
                            double ll_split_lambda,
                            bool ll_split_weight_penalty,
                            std::vector<size_t> ll_split_variables,
                            size_t ll_split_cutoff,
                            std::vector<double> overall_beta,
                            unsigned int mtry,
                            unsigned int num_trees,
                            unsigned int min_node_size,
                            double sample_fraction,
                            bool honesty,
                            double honesty_fraction,
                            bool honesty_prune_leaves,
                            size_t ci_group_size,
                            double alpha,
                            double imbalance_penalty,
                            std::vector<size_t> clusters,
                            unsigned int samples_per_cluster,
                            unsigned int num_threads,
                            unsigned int seed) {
  ForestTrainer trainer = ll_regression_trainer(ll_split_lambda, ll_split_weight_penalty, overall_beta,
                                               ll_split_cutoff, ll_split_variables);

  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_outcome_index(outcome_index);

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
                        honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(data, options);

  std::vector<Prediction> predictions;
  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List ll_regression_predict(Rcpp::List forest_object,
                                const Rcpp::NumericMatrix& train_matrix,
                                size_t outcome_index,
                                const Rcpp::NumericMatrix& test_matrix,
                                std::vector<double> ll_lambda,
                                bool ll_weight_penalty,
                                std::vector<size_t> linear_correction_variables,
                                unsigned int num_threads,
                                bool estimate_variance) {
  Data train_data = RcppUtilities::convert_data(train_matrix);
  train_data.set_outcome_index(outcome_index);
  Data data = RcppUtilities::convert_data(test_matrix);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = ll_regression_predictor(num_threads,
      ll_lambda, ll_weight_penalty, linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, train_data, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}

// [[Rcpp::export]]
Rcpp::List ll_regression_predict_oob(Rcpp::List forest_object,
                                    const Rcpp::NumericMatrix& train_matrix,
                                    size_t outcome_index,
                                    std::vector<double> ll_lambda,
                                    bool ll_weight_penalty,
                                    std::vector<size_t> linear_correction_variables,
                                    unsigned int num_threads,
                                    bool estimate_variance) {
  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_outcome_index(outcome_index);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = ll_regression_predictor(num_threads,
      ll_lambda, ll_weight_penalty, linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}
