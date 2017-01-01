#include <Rcpp.h>
#include <vector>
#include <sstream>
#include <map>

#include "globals.h"
#include "RcppUtilities.h"
#include "QuantileRelabelingStrategy.h"
#include "ProbabilitySplittingRule.h"
#include "QuantilePredictionStrategy.h"

Rcpp::List train(std::vector<double> &quantiles,
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

  RelabelingStrategy *relabeling_strategy = new QuantileRelabelingStrategy(&quantiles);
  SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, quantiles.size());
  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(&quantiles);

  std::unordered_map<std::string, size_t> observables = {{"outcome", outcome_index}};
  ForestModel *forest_model = new ForestModel(observables,
                                              relabeling_strategy,
                                              splitting_rule,
                                              prediction_strategy);

  RcppUtilities::initializeForestModel(forest_model, mtry, num_trees, num_threads,
      min_node_size, sample_with_replacement, sample_fraction);

  Forest *forest = forest_model->train(data);

  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);

  Rcpp::List result;
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest->get_trees()->size(), "num.trees");

  delete forest_model;
  delete forest;
  delete data;

  return result;
}


Rcpp::NumericMatrix predict(std::vector<double> &quantiles,
                            Rcpp::List forest,
                            Rcpp::NumericMatrix input_data,
                            uint outcome_index,
                            Rcpp::RawMatrix sparse_data,
                            std::vector <std::string> variable_names,
                            uint mtry,
                            uint num_trees,
                            bool verbose,
                            uint num_threads,
                            uint min_node_size,
                            bool sample_with_replacement,
                            bool keep_inbag,
                            double sample_fraction) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  RelabelingStrategy *relabeling_strategy = new QuantileRelabelingStrategy(&quantiles);
  SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, quantiles.size());
  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(&quantiles);

  std::unordered_map<std::string, size_t> observables = {{"outcome", outcome_index}};
  ForestModel *forest_model = new ForestModel(observables,
                                              relabeling_strategy,
                                              splitting_rule,
                                              prediction_strategy);

  RcppUtilities::initializeForestModel(forest_model, mtry, num_trees, num_threads,
                        min_node_size, sample_with_replacement, sample_fraction);

  Forest* deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  std::vector<std::vector<double>> predictions = forest_model->predict(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions, quantiles.size());

  delete forest_model;
  delete deserialized_forest;
  delete data;

  return result;
}

// [[Rcpp::export]]
Rcpp::List train_test(std::vector<double> &quantiles,
                      Rcpp::NumericMatrix input_data,
                      uint outcome_index) {
  return train(quantiles, input_data, outcome_index,
               Rcpp::RawMatrix(),
               std::vector<std::string>(),
               5, 5, false, 2, 5, true, false, 0.7);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix predict_test(Rcpp::List forest,
                                 std::vector<double> &quantiles,
                                 Rcpp::NumericMatrix input_data) {
  return predict(quantiles, forest, input_data,
                 -1,
                 Rcpp::RawMatrix(),
                 std::vector<std::string>(),
                 5, 5, false, 2, 5, true, false, 0.7);
}
