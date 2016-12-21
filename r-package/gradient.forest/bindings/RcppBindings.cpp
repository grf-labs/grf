#include <Rcpp.h>
#include <vector>
#include <sstream>
#include <map>

#include "globals.h"
#include "DataDouble.h"
#include "QuantileRelabelingStrategy.h"
#include "ProbabilitySplittingRule.h"
#include "QuantilePredictionStrategy.h"
#include "Forest.h"
#include "ForestModel.h"

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

  forest_model->initCpp(mtry, num_trees, &std::cout, 0, num_threads, load_forest_filename,
                        min_node_size, split_select_weights_file, always_split_variable_names,
                        sample_with_replacement,
                        memory_saving_splitting, case_weights_file, sample_fraction);
}

// [[Rcpp::export]]
Rcpp::List train(std::vector<double> &quantiles,
                 Rcpp::NumericMatrix input_data,
                 Rcpp::RawMatrix sparse_data,
                 uint outcome_index,
                 std::vector <std::string> variable_names,
                 uint mtry,
                 uint num_trees,
                 bool verbose,
                 uint num_threads,
                 uint min_node_size,
                 bool sample_with_replacement,
                 bool keep_inbag,
                 double sample_fraction) {

  Rcpp::List result;

  size_t num_rows = input_data.nrow();
  size_t num_cols = input_data.ncol();

  Data *data = new DataDouble(input_data.begin(), variable_names, num_rows, num_cols);

  if (sparse_data.nrow() > 1) {
    data->addSparseData(sparse_data.begin(), sparse_data.ncol());
  }

  std::vector<double> *initialized_quantiles = !quantiles.empty()
                                               ? new std::vector<double>(quantiles)
                                               : new std::vector<double>({0.25, 0.5, 0.75});
  RelabelingStrategy *relabeling_strategy = new QuantileRelabelingStrategy(&quantiles);
  SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, quantiles.size());
  PredictionStrategy *prediction_strategy = new QuantilePredictionStrategy(&quantiles);

  std::unordered_map<std::string, size_t> observables = {{"outcome", outcome_index}};
  ForestModel *forest_model = new ForestModel(observables,
                                              relabeling_strategy,
                                              splitting_rule,
                                              prediction_strategy);

  initializeForestModel(forest_model, mtry, num_trees, num_threads,
      min_node_size, sample_with_replacement, sample_fraction);


  Forest *forest = forest_model->train(data);

  if (keep_inbag) {
    result.push_back(forest->get_inbag_counts(), "inbag.counts");
  }

  result.push_back(forest->get_trees()->size(), "num.trees");

  delete forest_model;
  delete forest;
  delete data;

  return result;
}