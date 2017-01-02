#include <Rcpp.h>
#include <vector>
#include <sstream>
#include <map>

#include "globals.h"
#include "RcppUtilities.h"
#include "InstrumentalRelabelingStrategy.h"
#include "InstrumentalPredictionStrategy.h"
#include "RegressionSplittingRule.h"

// [[Rcpp::export]]
Rcpp::List instrumental_train(Rcpp::NumericMatrix input_data,
                              uint outcome_index,
                              uint treatment_index,
                              uint instrument_index,
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

  std::unordered_map<std::string, size_t> observables = {
      {Observations::OUTCOME, outcome_index},
      {Observations::TREATMENT, treatment_index},
      {Observations::INSTRUMENT, instrument_index}};

  RelabelingStrategy *relabeling_strategy = new InstrumentalRelabelingStrategy();
  SplittingRule *splitting_rule = new RegressionSplittingRule(data);
  PredictionStrategy *prediction_strategy = new InstrumentalPredictionStrategy();

  ForestTrainer *forest_trainer = new ForestTrainer(observables,
                                                    relabeling_strategy,
                                                    splitting_rule,
                                                    prediction_strategy);
  RcppUtilities::initialize_forest_trainer(forest_trainer, mtry, num_trees, num_threads,
      min_node_size, sample_with_replacement, sample_fraction, no_split_variables, seed);

  Forest *forest = forest_trainer->train(data);

  Rcpp::List result;
  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest->get_trees()->size(), "num.trees");

  delete forest_trainer;
  delete forest;
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

  PredictionStrategy *prediction_strategy = new InstrumentalPredictionStrategy();
  ForestPredictor* forest_predictor = new ForestPredictor(prediction_strategy);
  forest_predictor->init("", num_threads, &std::cout);

  Forest* deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  std::vector<std::vector<double>> predictions = forest_predictor->predict(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions, 1);

  delete forest_predictor;
  delete deserialized_forest;
  delete data;

  return result;
}
