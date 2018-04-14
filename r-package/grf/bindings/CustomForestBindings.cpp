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
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List custom_train(Rcpp::NumericMatrix input_data,
                        Eigen::SparseMatrix<double> sparse_input_data,
                        size_t outcome_index,
                        unsigned int mtry,
                        unsigned int num_trees,
                        unsigned int num_threads,
                        unsigned int min_node_size,
                        double sample_fraction,
                        unsigned int seed,
                        bool honesty,
                        unsigned int ci_group_size,
                        double alpha,
                        double imbalance_penalty,
                        std::vector<size_t> clusters,
                        unsigned int samples_per_cluster) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
      honesty, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

  ForestTrainer trainer = ForestTrainers::custom_trainer(outcome_index - 1);
  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix custom_predict(Rcpp::List forest_object,
                                   Rcpp::NumericMatrix input_data,
                                   Eigen::SparseMatrix<double> sparse_input_data,
                                   unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::custom_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix custom_predict_oob(Rcpp::List forest_object,
                                       Rcpp::NumericMatrix input_data,
                                       Eigen::SparseMatrix<double> sparse_input_data,
                                       unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::custom_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}
