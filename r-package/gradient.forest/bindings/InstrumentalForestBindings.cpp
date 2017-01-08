#include <Rcpp.h>
#include <vector>
#include <sstream>
#include <map>

#include "globals.h"
#include "RcppUtilities.h"

#include "ForestTrainers.h"
#include "ForestPredictors.h"

// [[Rcpp::export]]
Rcpp::List instrumental_train(Rcpp::NumericMatrix input_data,
                              size_t outcome_index,
                              size_t treatment_index,
                              size_t instrument_index,
                              Rcpp::RawMatrix sparse_data,
                              std::vector<std::string> variable_names,
                              uint mtry,
                              uint num_trees,
                              bool verbose,
                              uint num_threads,
                              uint min_node_size,
                              bool sample_with_replacement,
                              bool keep_inbag,
                              double sample_fraction,
                              std::vector<size_t> no_split_variables,
                              uint seed) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  ForestTrainer trainer = ForestTrainers::instrumental_trainer(data, outcome_index,
      treatment_index, instrument_index);
  RcppUtilities::initialize_trainer(trainer, mtry, num_trees, num_threads,
      min_node_size, sample_with_replacement, sample_fraction, no_split_variables, seed);

  Forest forest = trainer.train(data);

  Rcpp::List result;
  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest.get_trees().size(), "num.trees");

  delete data;

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix instrumental_predict(Rcpp::List forest,
                                         Rcpp::NumericMatrix input_data,
                                         Rcpp::RawMatrix sparse_data,
                                         std::vector <std::string> variable_names,
                                         uint num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads);
  std::vector<std::vector<double>> predictions = predictor.predict(deserialized_forest, data);

  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions, 1);

  delete data;

  return result;
}
