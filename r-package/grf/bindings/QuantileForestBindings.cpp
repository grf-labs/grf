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
Rcpp::List quantile_train(std::vector<double> quantiles,
                          bool regression_splitting,
                          const Rcpp::NumericMatrix& train_matrix,
                          size_t outcome_index,
                          unsigned int mtry,
                          unsigned int num_trees,
                          int min_node_size,
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
                          int num_threads,
                          unsigned int seed) {
  ForestTrainer trainer = regression_splitting
      ? regression_trainer()
      : quantile_trainer(quantiles);

  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_outcome_index(outcome_index);

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = quantile_predictor(num_threads, quantiles);
    predictions = predictor.predict_oob(forest, data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict(const Rcpp::List& forest_object,
                                     std::vector<double> quantiles,
                                     const Rcpp::NumericMatrix& train_matrix,
                                     size_t outcome_index,
                                     const Rcpp::NumericMatrix& test_matrix,
                                     unsigned int num_threads) {
  Data train_data = RcppUtilities::convert_data(train_matrix);
  Data data = RcppUtilities::convert_data(test_matrix);
  train_data.set_outcome_index(outcome_index);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict(forest, train_data, data, false);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict_oob(const Rcpp::List& forest_object,
                                         std::vector<double> quantiles,
                                         const Rcpp::NumericMatrix& train_matrix,
                                         size_t outcome_index,
                                         unsigned int num_threads) {
  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_outcome_index(outcome_index);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, false);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  return result;
}
