#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List quantile_train(std::vector<double> quantiles,
                          Rcpp::NumericMatrix input_data,
                          size_t outcome_index,
                          Rcpp::RawMatrix sparse_data,
                          std::vector<std::string> variable_names,
                          unsigned int mtry,
                          unsigned int num_trees,
                          bool verbose,
                          int num_threads,
                          int min_node_size,
                          bool sample_with_replacement,
                          bool keep_inbag,
                          double sample_fraction,
                          std::vector<size_t> no_split_variables,
                          unsigned int seed,
                          bool honesty,
                          unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  ForestTrainer trainer = ForestTrainers::quantile_trainer(data, outcome_index, quantiles);
  RcppUtilities::initialize_trainer(trainer, mtry, num_trees, num_threads, min_node_size,
      sample_with_replacement, sample_fraction, no_split_variables, seed, honesty, ci_group_size);
  Forest forest = trainer.train(data);

  Rcpp::List result;
  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest.get_trees().size(), "num.trees");

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict(Rcpp::List forest,
                                     std::vector<double> quantiles,
                                     Rcpp::NumericMatrix input_data,
                                     Rcpp::RawMatrix sparse_data,
                                     std::vector <std::string> variable_names,
                                     unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict_oob(Rcpp::List forest,
                                         std::vector<double> quantiles,
                                         Rcpp::NumericMatrix input_data,
                                         Rcpp::RawMatrix sparse_data,
                                         std::vector <std::string> variable_names,
                                         unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(num_threads, quantiles);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}
