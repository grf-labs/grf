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
                            std::vector<std::string> variable_names,
                            unsigned int mtry,
                            unsigned int num_trees,
                            bool verbose,
                            unsigned int num_threads,
                            unsigned int min_node_size,
                            bool sample_with_replacement,
                            bool keep_inbag,
                            double sample_fraction,
                            unsigned int seed,
                            bool honesty,
                            unsigned int ci_group_size,
                            double alpha,
                            double lambda,
                            bool downweight_penalty) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
                        honesty, sample_with_replacement, num_threads, seed);
  ForestTrainer trainer = lambda > 0
      ? ForestTrainers::regularized_regression_trainer(outcome_index - 1, lambda, downweight_penalty, options)
      : ForestTrainers::regression_trainer(outcome_index - 1, alpha, options);

  Forest forest = trainer.train(data);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List regression_predict(Rcpp::List forest_object,
                              Rcpp::NumericMatrix input_data,
                              Eigen::SparseMatrix<double> sparse_input_data,
                              std::vector<std::string> variable_names,
                              unsigned int num_threads,
                              unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);
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
                                  std::vector<std::string> variable_names,
                                  unsigned int num_threads,
                                  unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}
