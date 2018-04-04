#include <Rcpp.h>
#include <sstream>

#include "commons/DefaultData.h"
#include "commons/SparseData.h"
#include "forest/ForestOptions.h"
#include "RcppUtilities.h"
#include "serialization/ForestSerializer.h"

const std::string RcppUtilities::SERIALIZED_FOREST_KEY = "serialized.forest";

Rcpp::List RcppUtilities::create_forest_object(const Forest& forest,
                                               Data* data) {
  Rcpp::List result;
  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest.get_trees().size(), "num.trees");
  return result;
};


Rcpp::RawVector RcppUtilities::serialize_forest(const Forest& forest) {
  ForestSerializer forest_serializer;
  std::stringstream stream;
  forest_serializer.serialize(stream, forest);

  std::string contents = stream.str();

  Rcpp::RawVector result(contents.size());
  std::copy(contents.begin(), contents.end(), result.begin());
  return result;
}

Forest RcppUtilities::deserialize_forest(Rcpp::RawVector input) {
  ForestSerializer forest_serializer;

  std::string contents(input.begin(), input.end());

  std::stringstream stream(contents);
  return forest_serializer.deserialize(stream);
}

Data* RcppUtilities::convert_data(Rcpp::NumericMatrix input_data,
                                  Eigen::SparseMatrix<double>& sparse_input_data) {
  Data* data;

  if (input_data.nrow() > 0) {
    size_t num_rows = input_data.nrow();
    size_t num_cols = input_data.ncol();
    data = new DefaultData(input_data.begin(), num_rows, num_cols);
  } else {
    size_t num_rows = sparse_input_data.rows();
    size_t num_cols = sparse_input_data.cols();
    data = new SparseData(&sparse_input_data, num_rows, num_cols);
  }

  data->sort();
  return data;
}

Rcpp::List RcppUtilities::create_prediction_object(const std::vector<Prediction>& predictions) {
  Rcpp::List result;
  result.push_back(RcppUtilities::create_prediction_matrix(predictions), "predictions");
  result.push_back(RcppUtilities::create_variance_matrix(predictions), "variance.estimates");
  result.push_back(RcppUtilities::create_error_matrix(predictions), "debiased.error");
  return result;
};

Rcpp::NumericMatrix RcppUtilities::create_prediction_matrix(const std::vector<Prediction>& predictions) {
  if (predictions.empty()) {
    return Rcpp::NumericMatrix(0);
  }

  size_t prediction_length = predictions.at(0).size();
  Rcpp::NumericMatrix result(predictions.size(), prediction_length);

  for (size_t i = 0; i < predictions.size(); i++) {
    const std::vector<double>& prediction = predictions[i].get_predictions();
    for (size_t j = 0; j < prediction.size(); j++) {
      double value = prediction[j];
      result(i, j) = value;
    }
  }
  return result;
}

Rcpp::NumericMatrix RcppUtilities::create_variance_matrix(const std::vector<Prediction>& predictions) {
  if (predictions.empty()) {
    return Rcpp::NumericMatrix(0);
  }

  Prediction first_prediction = predictions.at(0);
  if (!first_prediction.contains_variance_estimates()) {
    return Rcpp::NumericMatrix(0);
  }

  size_t prediction_length = first_prediction.size();
  Rcpp::NumericMatrix result(predictions.size(), prediction_length);

  for (size_t i = 0; i < predictions.size(); i++) {
    const std::vector<double>& variance_estimate = predictions[i].get_variance_estimates();
    for (size_t j = 0; j < variance_estimate.size(); j++) {
      double value = variance_estimate[j];
      result(i, j) = value;
    }
  }
  return result;
}

Rcpp::NumericMatrix RcppUtilities::create_error_matrix(const std::vector<Prediction>& predictions) {
  if (predictions.empty()) {
    return Rcpp::NumericMatrix(0);
  }

  Prediction first_prediction = predictions.at(0);
  if (!first_prediction.contains_error_estimates()) {
    return Rcpp::NumericMatrix(0);
  }

  size_t prediction_length = first_prediction.size();
  Rcpp::NumericMatrix result(predictions.size(), prediction_length);

  for (size_t i = 0; i < predictions.size(); i++) {
    const std::vector<double>& error_estimate = predictions[i].get_error_estimates();
    for (size_t j = 0; j < error_estimate.size(); j++) {
      double value = error_estimate[j];
      result(i, j) = value;
    }
  }
  return result;
}
