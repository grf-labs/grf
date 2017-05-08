#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>


#include "commons/globals.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List instrumental_train(Rcpp::NumericMatrix input_data,
                              size_t outcome_index,
                              size_t treatment_index,
                              size_t instrument_index,
                              Rcpp::RawMatrix sparse_data,
                              std::vector<std::string> variable_names,
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
                              double split_regularization) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  ForestTrainer trainer = ForestTrainers::instrumental_trainer(data, outcome_index,
      treatment_index, instrument_index, split_regularization);
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
Rcpp::List instrumental_predict(Rcpp::List forest,
                                Rcpp::NumericMatrix input_data,
                                Rcpp::RawMatrix sparse_data,
                                std::vector <std::string> variable_names,
                                unsigned int num_threads,
                                unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, data);

  Rcpp::List result;
  result.push_back(RcppUtilities::create_prediction_matrix(predictions), "predictions");
  result.push_back(RcppUtilities::create_variance_matrix(predictions), "variance.estimates");

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List instrumental_predict_oob(Rcpp::List forest,
                                    Rcpp::NumericMatrix input_data,
                                    Rcpp::RawMatrix sparse_data,
                                    std::vector <std::string> variable_names,
                                    unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, 1);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data);

  Rcpp::List result;
  result.push_back(RcppUtilities::create_prediction_matrix(predictions), "predictions");
  result.push_back(RcppUtilities::create_variance_matrix(predictions), "variance.estimates");

  delete data;
  return result;
}
