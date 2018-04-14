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
                            unsigned int ci_group_size,
                            double alpha,
                            double imbalance_penalty,
                            std::vector<size_t> clusters,
                            unsigned int samples_per_cluster) {
  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index - 1);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
      honesty, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

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
                              unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List regression_predict_oob(Rcpp::List forest_object,
                                  Rcpp::NumericMatrix input_data,
                                  Eigen::SparseMatrix<double> sparse_input_data,
                                  unsigned int num_threads,
                                  unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
          forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

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
                                double lambda,
                                bool use_unweighted_penalty,
                                unsigned int num_threads) {
  Data *test_data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Data *original_data = RcppUtilities::convert_data(training_data, sparse_training_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::local_linear_predictor(num_threads, original_data, test_data,
                                                                       lambda, use_unweighted_penalty);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, test_data);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete original_data;
  delete test_data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List local_linear_predict_oob(Rcpp::List forest,
                                    Rcpp::NumericMatrix input_data,
                                    Eigen::SparseMatrix<double> sparse_input_data,
                                    double lambda,
                                    bool use_unweighted_penalty,
                                    unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);

  Forest deserialized_forest = RcppUtilities::deserialize_forest(forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::local_linear_predictor(num_threads, data, data,
                                                                       lambda, use_unweighted_penalty);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}

