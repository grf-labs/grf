#include <Rcpp.h>
#include <sstream>

#include "ForestSerializer.h"
#include "RcppUtilities.h"

const std::string RcppUtilities::SERIALIZED_FOREST_KEY = "serialized.forest";

void RcppUtilities::initialize_forest_trainer(ForestTrainer &forest_trainer,
                                              uint mtry,
                                              uint num_trees,
                                              uint num_threads,
                                              uint min_node_size,
                                              bool sample_with_replacement,
                                              double sample_fraction,
                                              std::vector<size_t> no_split_variables,
                                              uint seed) {
  std::string load_forest_filename = "";
  std::string split_select_weights_file = "";
  std::vector<std::string> always_split_variable_names;
  bool memory_saving_splitting = false;
  std::string case_weights_file = "";

  forest_trainer.init(mtry, num_trees, &std::cout, seed, num_threads, load_forest_filename,
                       min_node_size, no_split_variables, split_select_weights_file,
                       always_split_variable_names, sample_with_replacement,
                       memory_saving_splitting, case_weights_file, sample_fraction);
}

Rcpp::RawVector RcppUtilities::serialize_forest(Forest* forest) {
  ForestSerializer forest_serializer;
  std::stringstream stream;
  forest_serializer.serialize(stream, forest);

  std::string contents = stream.str();

  Rcpp::RawVector result(contents.size());
  std::copy(contents.begin(), contents.end(), result.begin());
  return result;
}

Forest* RcppUtilities::deserialize_forest(Rcpp::RawVector input) {
  ForestSerializer forest_serializer;

  std::string contents(input.begin(), input.end());

  std::stringstream stream(contents);
  return forest_serializer.deserialize(stream);
}

Data* RcppUtilities::convert_data(Rcpp::NumericMatrix input_data,
                                  Rcpp::RawMatrix sparse_data,
                                  std::vector<std::string> variable_names) {
  size_t num_rows = input_data.nrow();
  size_t num_cols = input_data.ncol();

  Data* data = new Data(input_data.begin(), variable_names, num_rows, num_cols);
  data->sort();

  if (sparse_data.nrow() > 1) {
    data->addSparseData(sparse_data.begin(), sparse_data.ncol());
  }

  return data;
}

Rcpp::NumericMatrix RcppUtilities::create_prediction_matrix(std::vector<std::vector<double>> predictions,
                                                            size_t prediction_length) {
  Rcpp::NumericMatrix result(predictions.size(), prediction_length);
  for (int i = 0; i < predictions.size(); i++) {
    std::vector<double> prediction = predictions[i];
    if (prediction.size() != prediction_length) {
      throw std::runtime_error("Prediction " + std::to_string(i) + " did not have the expected length.");
    }

    std::copy(prediction.begin(),
              prediction.end(),
              result.begin() + i * prediction_length);
  }
  return result;
}