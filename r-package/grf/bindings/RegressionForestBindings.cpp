#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "tuning/ParameterTuner.h"
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
                            bool downweight_penalty,
                            bool tune_parameters) {
  ForestTrainer trainer = lambda > 0
      ? ForestTrainers::regularized_regression_trainer(outcome_index - 1, lambda, downweight_penalty)
      : ForestTrainers::regression_trainer(outcome_index - 1, alpha);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
                        honesty, sample_with_replacement, num_threads, seed);

  if (tune_parameters) {
    ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads, ci_group_size);
    ParameterTuner tuner(trainer, predictor, outcome_index);
    uint tuned_min_node_size = tuner.tune_min_node_size(data, options);
    options.set_min_node_size(tuned_min_node_size);
  }

  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  result.push_back(options.get_min_node_size(), "tuned.min.node.size");

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
