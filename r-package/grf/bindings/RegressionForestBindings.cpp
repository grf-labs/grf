#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List regression_train(Rcpp::NumericMatrix input_data,
                            size_t outcome_index,
                            std::vector <std::string> variable_names,
                            unsigned int mtry,
                            unsigned int num_trees,
                            bool verbose,
                            unsigned int num_threads,
                            unsigned int min_node_size,
                            bool sample_with_replacement,
                            bool keep_inbag,
                            double sample_fraction,
                            std::vector<size_t> no_split_variables,
                            unsigned int seed,
                            bool honesty,
                            unsigned int ci_group_size,
                            double alpha,
                            double lambda,
                            bool downweight_penalty) {
  Data* data = RcppUtilities::convert_data(input_data, variable_names);

  ForestTrainer trainer = lambda > 0
      ? ForestTrainers::regularized_regression_trainer(data, outcome_index - 1, lambda, downweight_penalty)
      : ForestTrainers::regression_trainer(data, outcome_index - 1, alpha);

  RcppUtilities::initialize_trainer(trainer, mtry, num_trees, num_threads, min_node_size,
      sample_with_replacement, sample_fraction, no_split_variables, seed, honesty, ci_group_size);
  Forest forest = trainer.train(data);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List regression_predict(Rcpp::List forest_object,
                              Rcpp::NumericMatrix input_data,
                              std::vector<std::string> variable_names,
                              unsigned int num_threads,
                              unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, variable_names);
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
                                  std::vector<std::string> variable_names,
                                  unsigned int num_threads,
                                  unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, variable_names);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::regression_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}
