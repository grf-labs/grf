#include <Rcpp.h>
#include <sstream>

#include "commons/DefaultData.h"
#include "commons/SparseData.h"
#include "forest/ForestOptions.h"
#include "RcppUtilities.h"

Forest RcppUtilities::deserialize_forest(Rcpp::List forest_object) {
  size_t ci_group_size = forest_object["_ci_group_size"];
  size_t num_variables = forest_object["_num_variables"];

  size_t num_trees = forest_object["num.trees"];
  std::vector<std::shared_ptr<Tree>> trees(num_trees);

  std::vector<size_t> root_nodes = forest_object["_root_nodes"];
  std::vector<std::vector<std::vector<size_t>>> child_nodes = forest_object["_child_nodes"];
  std::vector<std::vector<std::vector<size_t>>> leaf_samples = forest_object["_leaf_samples"];

  std::vector<std::vector<size_t>> split_vars = forest_object["_split_vars"];
  std::vector<std::vector<double>> split_values = forest_object["_split_values"];
  std::vector<std::vector<size_t>> drawn_samples = forest_object["_drawn_samples"];

  std::vector<std::vector<std::vector<double>>> prediction_values = forest_object["_pv_values"];
  size_t num_types = forest_object["_pv_num_types"];

  for (size_t t = 0; t < num_trees; t++) {
    trees[t] = std::shared_ptr<Tree>(
        new Tree(root_nodes.at(t),
                 child_nodes.at(t),
                 leaf_samples.at(t),
                 split_vars.at(t),
                 split_values.at(t),
                 drawn_samples.at(t),
                 PredictionValues(prediction_values.at(t), num_types)));
  }

  return Forest(trees, num_variables, ci_group_size);
}

Rcpp::List RcppUtilities::serialize_forest(const Forest& forest) {
  Rcpp::List result;

  result.push_back(forest.get_ci_group_size(), "_ci_group_size");
  result.push_back(forest.get_num_variables(), "_num_variables");

  size_t num_trees = forest.get_trees().size();
  result.push_back(num_trees, "num.trees");

  std::vector<size_t> root_nodes(num_trees);
  std::vector<std::vector<std::vector<size_t>>> child_nodes(num_trees);
  std::vector<std::vector<std::vector<size_t>>> leaf_samples(num_trees);
  std::vector<std::vector<size_t>> split_vars(num_trees);
  std::vector<std::vector<double>> split_values(num_trees);
  std::vector<std::vector<size_t>> drawn_samples(num_trees);
  std::vector<std::vector<std::vector<double>>> prediction_values(num_trees);
  size_t num_types = 0;

  for (size_t t = 0; t < num_trees; t++) {
    std::shared_ptr<Tree> tree = forest.get_trees().at(t);
    root_nodes[t] = tree->get_root_node();
    child_nodes[t] = tree->get_child_nodes();
    leaf_samples[t] = tree->get_leaf_samples();
    split_vars[t] = tree->get_split_vars();
    split_values[t] = tree->get_split_values();
    drawn_samples[t] = tree->get_drawn_samples();

    prediction_values[t] = tree->get_prediction_values().get_all_values();
    if (t == 0) {
      num_types = tree->get_prediction_values().get_num_types();
    }
  }

  result.push_back(root_nodes, "_root_nodes");
  result.push_back(child_nodes, "_child_nodes");
  result.push_back(leaf_samples, "_leaf_samples");
  result.push_back(split_vars, "_split_vars");
  result.push_back(split_values, "_split_values");
  result.push_back(drawn_samples, "_drawn_samples");
  result.push_back(prediction_values, "_pv_values");
  result.push_back(num_types, "_pv_num_types");
  return result;
};

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
  return data;
}

Rcpp::List RcppUtilities::create_prediction_object(const std::vector<Prediction>& predictions) {
  Rcpp::List result;
  result.push_back(RcppUtilities::create_prediction_matrix(predictions), "predictions");
  result.push_back(RcppUtilities::create_variance_matrix(predictions), "variance.estimates");
  result.push_back(RcppUtilities::create_error_matrix(predictions), "debiased.error");
  result.push_back(RcppUtilities::create_excess_error_matrix(predictions), "excess.error");
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


Rcpp::NumericMatrix RcppUtilities::create_excess_error_matrix(const std::vector<Prediction>& predictions) {
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
    const std::vector<double>& error_estimate = predictions[i].get_excess_error_estimates();
    for (size_t j = 0; j < error_estimate.size(); j++) {
      double value = error_estimate[j];
      result(i, j) = value;
    }
  }
  return result;
}

