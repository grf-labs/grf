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
Rcpp::List regression_train(Rcpp::NumericMatrix input_data,
                            Eigen::SparseMatrix<double> sparse_input_data,
                            size_t outcome_index,
                            unsigned int mtry,
                            unsigned int num_trees,
                            unsigned int num_threads,
                            unsigned int min_node_size,
                            double sample_fraction,
                            unsigned int seed,
                            bool honesty,
                            double honesty_fraction,
                            size_t ci_group_size,
                            double alpha,
                            double imbalance_penalty,
                            std::vector<size_t> clusters,
                            unsigned int samples_per_cluster) {
  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index - 1);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  result.push_back(options.get_tree_options().get_min_node_size(), "min.node.size");

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List regression_predict(Rcpp::List forest_object,
                              Rcpp::NumericMatrix input_data,
                              Eigen::SparseMatrix<double> sparse_input_data,
                              unsigned int num_threads,
                              bool estimate_variance) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, data, estimate_variance);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List regression_predict_oob(Rcpp::List forest_object,
                                  Rcpp::NumericMatrix input_data,
                                  Eigen::SparseMatrix<double> sparse_input_data,
                                  unsigned int num_threads,
                                  bool estimate_variance) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, estimate_variance);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List local_linear_predict(Rcpp::List forest,
                                Rcpp::NumericMatrix input_data,
                                Rcpp::NumericMatrix training_data,
                                Eigen::SparseMatrix<double> sparse_input_data,
                                Eigen::SparseMatrix<double> sparse_training_data,
                                std::vector<double> lambdas,
                                bool weight_penalty,
                                std::vector<size_t> linear_correction_variables,
                                unsigned int num_threads,
                                bool estimate_variance) {
  Data *test_data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Data *original_data = RcppUtilities::convert_data(training_data, sparse_training_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::local_linear_predictor(num_threads, original_data, test_data,
                                                                       lambdas, weight_penalty,
                                                                       linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, test_data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete original_data;
  delete test_data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List local_linear_predict_oob(Rcpp::List forest,
                                    Rcpp::NumericMatrix input_data,
                                    Eigen::SparseMatrix<double> sparse_input_data,
                                    std::vector<double> lambdas,
                                    bool weight_penalty,
                                    std::vector<size_t> linear_correction_variables,
                                    unsigned int num_threads,
                                    bool estimate_variance) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::local_linear_predictor(num_threads, data, data,
                                                                       lambdas, weight_penalty,
                                                                       linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}

