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
Rcpp::List quantile_train(std::vector<double> quantiles,
                          bool regression_splits,
                          Rcpp::NumericMatrix input_data,
                          Eigen::SparseMatrix<double> sparse_input_data,
                          size_t outcome_index,
                          unsigned int mtry,
                          unsigned int num_trees,
                          int num_threads,
                          int min_node_size,
                          double sample_fraction,
                          unsigned int seed,
                          bool honesty,
                          double honesty_fraction,
                          unsigned int ci_group_size,
                          double alpha,
                          double imbalance_penalty,
                          std::vector<size_t> clusters,
                          unsigned int samples_per_cluster) {
  ForestTrainer trainer = regression_splits
      ? ForestTrainers::regression_trainer(outcome_index - 1)
      : ForestTrainers::quantile_trainer(outcome_index - 1, quantiles);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict(Rcpp::List forest_object,
                                     std::vector<double> quantiles,
                                     Rcpp::NumericMatrix input_data,
                                     Eigen::SparseMatrix<double> sparse_input_data,
                                     unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict(forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict_oob(Rcpp::List forest_object,
                                         std::vector<double> quantiles,
                                         Rcpp::NumericMatrix input_data,
                                         Eigen::SparseMatrix<double> sparse_input_data,
                                         unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}
