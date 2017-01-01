#include <Rcpp.h>
#include <vector>
#include <sstream>
#include <map>

#include "globals.h"
#include "RcppUtilities.h"
#include "QuantileRelabelingStrategy.h"
#include "ProbabilitySplittingRule.h"
#include "QuantilePredictionStrategy.h"

// [[Rcpp::export]]
Rcpp::List quantile_train(std::vector<double> &quantiles,
                          Rcpp::NumericMatrix input_data,
                          uint outcome_index,
                          Rcpp::RawMatrix sparse_data,
                          std::vector<std::string> variable_names,
                          uint mtry,
                          uint num_trees,
                          bool verbose,
                          uint num_threads,
                          uint min_node_size,
                          bool sample_with_replacement,
                          bool keep_inbag,
                          double sample_fraction) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  std::unordered_map<std::string, size_t> observables = {{"outcome", outcome_index}};
  RelabelingStrategy *relabeling_strategy = new QuantileRelabelingStrategy(&quantiles);
  SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, quantiles.size());
  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(&quantiles);

  ForestTrainer *forest_trainer = new ForestTrainer(observables,
                                                    relabeling_strategy,
                                                    splitting_rule,
                                                    prediction_strategy);
  RcppUtilities::initialize_forest_trainer(forest_trainer, mtry, num_trees, num_threads,
      min_node_size, sample_with_replacement, sample_fraction);

  Forest *forest = forest_trainer->train(data);

  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);

  Rcpp::List result;
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest->get_trees()->size(), "num.trees");

  delete forest_trainer;
  delete forest;
  delete data;

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix quantile_predict(Rcpp::List forest,
                                     std::vector<double> &quantiles,
                                     Rcpp::NumericMatrix input_data,
                                     Rcpp::RawMatrix sparse_data,
                                     std::vector <std::string> variable_names,
                                     uint num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(&quantiles);

  ForestPredictor* forest_predictor = new ForestPredictor(prediction_strategy);
  forest_predictor->init("", num_threads, &std::cout);

  Forest* deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  std::vector<std::vector<double>> predictions = forest_predictor->predict(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions, quantiles.size());

  delete forest_predictor;
  delete deserialized_forest;
  delete data;

  return result;
}
