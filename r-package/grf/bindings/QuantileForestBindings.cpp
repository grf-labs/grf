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
                          Rcpp::NumericMatrix train_matrix,
                          Eigen::SparseMatrix<double> sparse_train_matrix,
                          size_t outcome_index,
                          unsigned int mtry,
                          unsigned int num_trees,
                          int min_node_size,
                          double sample_fraction,
                          bool honesty,
                          double honesty_fraction,
                          size_t ci_group_size,
                          double alpha,
                          double imbalance_penalty,
                          std::vector<size_t> clusters,
                          unsigned int samples_per_cluster,
                          int num_threads,
                          unsigned int seed) {
  ForestTrainer trainer = regression_splits
      ? ForestTrainers::regression_trainer()
      : ForestTrainers::quantile_trainer(quantiles);

  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->sort();

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(data, options);

  delete data;
  return RcppUtilities::create_forest_object(forest, std::vector<Prediction>());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict(Rcpp::List forest_object,
                                     std::vector<double> quantiles,
                                     Rcpp::NumericMatrix train_matrix,
                                     Eigen::SparseMatrix<double> sparse_train_matrix,
                                     size_t outcome_index,
                                     Rcpp::NumericMatrix test_matrix,
                                     Eigen::SparseMatrix<double> sparse_test_matrix,
                                     unsigned int num_threads) {
  Data* train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  Data* data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);
  train_data->set_outcome_index(outcome_index - 1);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict(forest, train_data, data, false);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete train_data;
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict_oob(Rcpp::List forest_object,
                                         std::vector<double> quantiles,
                                         Rcpp::NumericMatrix train_matrix,
                                         Eigen::SparseMatrix<double> sparse_train_matrix,
                                         size_t outcome_index,
                                         unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, false);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}
