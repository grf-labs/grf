#include <Rcpp.h>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List probability_train(Rcpp::NumericMatrix& train_matrix,
                             Eigen::SparseMatrix<double>& sparse_train_matrix,
                             size_t outcome_index,
                             size_t sample_weight_index,
                             bool use_sample_weights,
                             size_t num_classes,
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
                             const std::vector<size_t>& clusters,
                             unsigned int samples_per_cluster,
                             bool compute_oob_predictions,
                             int num_threads,
                             unsigned int seed) {
  ForestTrainer trainer = probability_trainer(num_classes);

  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index);
  if (use_sample_weights) {
      data->set_weight_index(sample_weight_index);
  }

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(*data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = probability_predictor(num_threads, num_classes);
    predictions = predictor.predict_oob(forest, *data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List probability_predict(Rcpp::List& forest_object,
                               Rcpp::NumericMatrix& train_matrix,
                               Eigen::SparseMatrix<double>& sparse_train_matrix,
                               size_t outcome_index,
                               size_t num_classes,
                               Rcpp::NumericMatrix& test_matrix,
                               Eigen::SparseMatrix<double>& sparse_test_matrix,
                               unsigned int num_threads,
                               bool estimate_variance) {
  std::unique_ptr<Data> train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  std::unique_ptr<Data> data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);
  train_data->set_outcome_index(outcome_index);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = probability_predictor(num_threads, num_classes);
  std::vector<Prediction> predictions = predictor.predict(forest, *train_data, *data, estimate_variance);

  return RcppUtilities::create_prediction_object(predictions);
}

// [[Rcpp::export]]
Rcpp::List probability_predict_oob(Rcpp::List& forest_object,
                                   Rcpp::NumericMatrix& train_matrix,
                                   Eigen::SparseMatrix<double>& sparse_train_matrix,
                                   size_t outcome_index,
                                   size_t num_classes,
                                   unsigned int num_threads,
                                   bool estimate_variance) {
  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = probability_predictor(num_threads, num_classes);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, *data, estimate_variance);

  return RcppUtilities::create_prediction_object(predictions);
}
