#include <Rcpp.h>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List multi_causal_train(Rcpp::NumericMatrix train_matrix,
                              Eigen::SparseMatrix<double> sparse_train_matrix,
                              size_t outcome_index,
                              const std::vector<size_t>& treatment_index,
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
  size_t num_treatments = treatment_index.size();
  ForestTrainer trainer = multi_causal_trainer(num_treatments);

  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index);
  data->set_treatment_index(treatment_index);
  if (use_sample_weights) {
    data->set_weight_index(sample_weight_index);
  }

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(*data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = multi_causal_predictor(num_threads, num_treatments);
    predictions = predictor.predict_oob(forest, *data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List multi_causal_predict(Rcpp::List forest_object,
                                Rcpp::NumericMatrix train_matrix,
                                Eigen::SparseMatrix<double> sparse_train_matrix,
                                Rcpp::NumericMatrix test_matrix,
                                Eigen::SparseMatrix<double> sparse_test_matrix,
                                size_t num_treatments,
                                unsigned int num_threads,
                                bool estimate_variance) {
  std::unique_ptr<Data> train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  std::unique_ptr<Data> data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = multi_causal_predictor(num_threads, num_treatments);
  std::vector<Prediction> predictions = predictor.predict(forest, *train_data, *data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}

// [[Rcpp::export]]
Rcpp::List multi_causal_predict_oob(Rcpp::List forest_object,
                                    Rcpp::NumericMatrix train_matrix,
                                    Eigen::SparseMatrix<double> sparse_train_matrix,
                                    size_t num_treatments,
                                    unsigned int num_threads,
                                    bool estimate_variance) {
  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = multi_causal_predictor(num_threads, num_treatments);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, *data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}
