#include <Rcpp.h>
#include <vector>
#include <sstream>
#include <map>

#include "globals.h"
#include "QuantileRelabelingStrategy.h"
#include "ProbabilitySplittingRule.h"
#include "QuantilePredictionStrategy.h"
#include "Forest.h"
#include "ForestModel.h"
#include "ForestSerializer.h"

void initializeForestModel(ForestModel *forest_model,
                           uint mtry,
                           uint num_trees,
                           uint num_threads,
                           uint min_node_size,
                           bool sample_with_replacement,
                           double sample_fraction) {
  std::string load_forest_filename = "";
  std::string split_select_weights_file = "";
  std::vector<std::string> always_split_variable_names;
  bool memory_saving_splitting = false;
  std::string case_weights_file = "";

  forest_model->init(mtry, num_trees, &std::cout, 0, num_threads, load_forest_filename,
                        min_node_size, split_select_weights_file, always_split_variable_names,
                        sample_with_replacement,
                        memory_saving_splitting, case_weights_file, sample_fraction);
}

Rcpp::RawVector serialize_forest(Forest* forest) {
  ForestSerializer forest_serializer;
  std::stringstream stream;
  forest_serializer.serialize(stream, forest);

  std::string contents = stream.str();

  Rcpp::RawVector result(contents.size());
  std::copy(contents.begin(), contents.end(), result.begin());
  return result;
}

Forest* deserialize_forest(Rcpp::RawVector input) {
  ForestSerializer forest_serializer;

  std::string contents(input.begin(), input.end());

  std::stringstream stream(contents);
  return forest_serializer.deserialize(stream);
}

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
  size_t num_rows = input_data.nrow();
  size_t num_cols = input_data.ncol();

  Data *data = new Data(input_data.begin(), variable_names, num_rows, num_cols);
  data->sort();

  if (sparse_data.nrow() > 1) {
    data->addSparseData(sparse_data.begin(), sparse_data.ncol());
  }

  std::vector<double> *initialized_quantiles = !quantiles.empty()
                                               ? new std::vector<double>(quantiles)
                                               : new std::vector<double>({0.25, 0.5, 0.75});
  RelabelingStrategy *relabeling_strategy = new QuantileRelabelingStrategy(&quantiles);
  SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, initialized_quantiles->size());
  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(initialized_quantiles);

  std::unordered_map<std::string, size_t> observables = {{"outcome", outcome_index}};
  ForestModel *forest_model = new ForestModel(observables,
                                              relabeling_strategy,
                                              splitting_rule,
                                              prediction_strategy);

  initializeForestModel(forest_model, mtry, num_trees, num_threads,
      min_node_size, sample_with_replacement, sample_fraction);

  Forest *forest = forest_model->train(data);

  Rcpp::RawVector serialized_forest = serialize_forest(forest);

  Rcpp::List result;
  result.push_back(serialized_forest, "serialized.forest");
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
  size_t num_rows = input_data.nrow();
  size_t num_cols = input_data.ncol();

  Data *data = new Data(input_data.begin(), variable_names, num_rows, num_cols);
  data->sort();

  if (sparse_data.nrow() > 1) {
    data->addSparseData(sparse_data.begin(), sparse_data.ncol());
  }

  std::vector<double> *initialized_quantiles = !quantiles.empty()
                                               ? new std::vector<double>(quantiles)
                                               : new std::vector<double>({0.25, 0.5, 0.75});
  RelabelingStrategy *relabeling_strategy = new QuantileRelabelingStrategy(&quantiles);
  SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, initialized_quantiles->size());
  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(initialized_quantiles);

  std::unordered_map<std::string, size_t> observables = {{"outcome", outcome_index}};
  ForestModel *forest_model = new ForestModel(observables,
                                              relabeling_strategy,
                                              splitting_rule,
                                              prediction_strategy);

  initializeForestModel(forest_model, mtry, num_trees, num_threads,
                        min_node_size, sample_with_replacement, sample_fraction);

  Forest* deserialized_forest = deserialize_forest(forest["serialized.forest"]);

  size_t num_quantiles = initialized_quantiles->size();
  std::vector<std::vector<double>> predictions = forest_model->predict(deserialized_forest, data);

  Rcpp::NumericMatrix result(predictions.size(), num_quantiles);
  for (int i = 0; i < predictions.size(); i++) {
    std::vector<double> prediction = predictions[i];
    if (prediction.size() != num_quantiles) {
      throw std::runtime_error("Prediction " + std::to_string(i) +
          " did not include a value for each quantile.");
    }

    std::copy(prediction.begin(),
              prediction.end(),
              result.begin() + i * num_quantiles);
  }


  // TODO: copy predictions

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
                                 Rcpp::NumericMatrix input_data,
                                 uint outcome_index) {
  return predict(quantiles, forest, input_data, outcome_index,
                 Rcpp::RawMatrix(),
                 std::vector<std::string>(),
                 5, 5, false, 2, 5, true, false, 0.7);
}
