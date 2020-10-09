#include <Rcpp.h>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List causal_survival_train(Rcpp::NumericMatrix& train_matrix,
                                 Eigen::SparseMatrix<double>& sparse_train_matrix,
                                 size_t causal_survival_numerator_index,
                                 size_t causal_survival_denominator_index,
                                 size_t treatment_index,
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
                                 size_t ci_group_size,
                                 double alpha,
                                 double imbalance_penalty,
                                 bool stabilize_splits,
                                 const std::vector<size_t>& clusters,
                                 unsigned int samples_per_cluster,
                                 bool compute_oob_predictions,
                                 unsigned int num_threads,
                                 unsigned int seed) {
  ForestTrainer trainer = causal_survival_trainer(stabilize_splits);

  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_causal_survival_numerator_index(causal_survival_numerator_index);
  data->set_causal_survival_denominator_index(causal_survival_denominator_index);
  data->set_treatment_index(treatment_index);
  data->set_instrument_index(treatment_index);
  data->set_censor_index(censor_index);
  if (use_sample_weights) {
    data->set_weight_index(sample_weight_index);
  }

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(*data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = causal_survival_predictor(num_threads);
    predictions = predictor.predict_oob(forest, *data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List causal_survival_predict(Rcpp::List& forest_object,
                                   Rcpp::NumericMatrix& train_matrix,
                                   Eigen::SparseMatrix<double>& sparse_train_matrix,
                                   Rcpp::NumericMatrix& test_matrix,
                                   Eigen::SparseMatrix<double>& sparse_test_matrix,
                                   unsigned int num_threads,
                                   bool estimate_variance) {
  std::unique_ptr<Data> train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  std::unique_ptr<Data> data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = causal_survival_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, *train_data, *data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}

// [[Rcpp::export]]
Rcpp::List causal_survival_predict_oob(Rcpp::List& forest_object,
                                       Rcpp::NumericMatrix& train_matrix,
                                       Eigen::SparseMatrix<double>& sparse_train_matrix,
                                       unsigned int num_threads,
                                       bool estimate_variance) {
  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = causal_survival_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, *data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}
