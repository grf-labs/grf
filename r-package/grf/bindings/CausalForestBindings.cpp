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
Rcpp::List causal_train(Rcpp::NumericMatrix input_data,
                        Eigen::SparseMatrix<double> sparse_input_data,
                        size_t outcome_index,
                        size_t treatment_index,
                        size_t instrument_index,
                        unsigned int mtry,
                        unsigned int num_trees,
                        unsigned int num_threads,
                        unsigned int min_node_size,
                        double sample_fraction,
                        unsigned int seed,
                        bool honesty,
                        unsigned int ci_group_size,
                        double reduced_form_weight,
                        double alpha,
                        bool imbalance_penalty,
                        bool stabilize_splits,
                        std::vector<size_t> clusters,
                        unsigned int samples_per_cluster) {
  ForestTrainer trainer = ForestTrainers::causal_trainer(
      outcome_index - 1,
      treatment_index - 1,
      instrument_index - 1,
      reduced_form_weight,
      stabilize_splits);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
      honesty, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}


// [[Rcpp::export]]
Rcpp::List causal_predict(Rcpp::List forest_object,
                          Rcpp::NumericMatrix input_data,
                          Eigen::SparseMatrix<double> sparse_input_data,
                          unsigned int num_threads,
                          unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List causal_predict_oob(Rcpp::List forest_object,
                              Rcpp::NumericMatrix input_data,
                              Eigen::SparseMatrix<double> sparse_input_data,
                              unsigned int num_threads,
                              unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List causal_predict_linear(Rcpp::List forest,
                                 Rcpp::NumericMatrix input_data,
                                 Rcpp::NumericMatrix training_data,
                                 Eigen::SparseMatrix<double> sparse_input_data,
                                 Eigen::SparseMatrix<double> sparse_training_data,
                                 std::vector<double> lambdas,
                                 bool use_unweighted_penalty,
                                 std::vector<size_t> linear_correction_variables,
                                 unsigned int num_threads) {
  Data *test_data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Data *original_data = RcppUtilities::convert_data(training_data, sparse_training_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::ll_causal_predictor(num_threads, original_data, test_data,
                                                                 lambdas, use_unweighted_penalty,
                                                                 linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, test_data);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete original_data;
  delete test_data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List causal_predict_oob_linear(Rcpp::List forest,
                                     Rcpp::NumericMatrix input_data,
                                     Eigen::SparseMatrix<double> sparse_input_data,
                                     std::vector<double> lambdas,
                                     bool use_unweighted_penalty,
                                     std::vector<size_t> linear_correction_variables,
                                     unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::ll_causal_predictor(num_threads, data, data,
                                                                 lambdas, use_unweighted_penalty,
                                                                 linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}
