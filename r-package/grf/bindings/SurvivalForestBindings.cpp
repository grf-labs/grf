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
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List survival_train(Rcpp::NumericMatrix train_matrix,
                          Eigen::SparseMatrix<double> sparse_train_matrix,
                          size_t outcome_index,
                          size_t censor_index,
                          size_t sample_weight_index,
                          bool use_sample_weights,
                          unsigned int mtry,
                          unsigned int num_trees,
                          unsigned int min_node_size,
                          double sample_fraction,
                          bool honesty,
                          double honesty_fraction,
                          bool honesty_prune_leaves,
                          double alpha,
                          size_t num_failures,
                          std::vector<size_t> clusters,
                          unsigned int samples_per_cluster,
                          bool compute_oob_predictions,
                          unsigned int num_threads,
                          unsigned int seed) {
  ForestTrainer trainer = survival_trainer();

  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_censor_index(censor_index - 1);
  if (use_sample_weights) {
      data->set_weight_index(sample_weight_index - 1);
  }

  size_t ci_group_size = 1;
  size_t imbalance_penalty = 0;
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(*data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = survival_predictor(num_threads, num_failures);
    predictions = predictor.predict_oob(forest, *data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List survival_predict(Rcpp::List forest_object,
                            Rcpp::NumericMatrix train_matrix,
                            Eigen::SparseMatrix<double> sparse_train_matrix,
                            size_t outcome_index,
                            size_t censor_index,
                            size_t sample_weight_index,
                            bool use_sample_weights,
                            Rcpp::NumericMatrix test_matrix,
                            Eigen::SparseMatrix<double> sparse_test_matrix,
                            unsigned int num_threads,
                            size_t num_failures) {
  std::unique_ptr<Data> train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  train_data->set_outcome_index(outcome_index - 1);
  train_data->set_censor_index(censor_index - 1);
  if (use_sample_weights) {
      train_data->set_weight_index(sample_weight_index - 1);
  }

  std::unique_ptr<Data> data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);
  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  bool estimate_variance = false;
  ForestPredictor predictor = survival_predictor(num_threads, num_failures);
  std::vector<Prediction> predictions = predictor.predict(forest, *train_data, *data, estimate_variance);

  return RcppUtilities::create_prediction_object(predictions);
}

// [[Rcpp::export]]
Rcpp::List survival_predict_oob(Rcpp::List forest_object,
                                Rcpp::NumericMatrix train_matrix,
                                Eigen::SparseMatrix<double> sparse_train_matrix,
                                size_t outcome_index,
                                size_t censor_index,
                                size_t sample_weight_index,
                                bool use_sample_weights,
                                unsigned int num_threads,
                                size_t num_failures) {
  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_censor_index(censor_index - 1);
  if (use_sample_weights) {
      data->set_weight_index(sample_weight_index - 1);
  }

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  bool estimate_variance = false;
  ForestPredictor predictor = survival_predictor(num_threads, num_failures);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, *data, estimate_variance);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  return result;
}
