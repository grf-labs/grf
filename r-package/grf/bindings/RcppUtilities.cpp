/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <Rcpp.h>

#include "commons/Data.h"
#include "forest/ForestOptions.h"
#include "RcppUtilities.h"

using namespace grf;

Rcpp::List RcppUtilities::create_forest_object(Forest& forest,
                                               const std::vector<Prediction>& predictions) {
  Rcpp::List result = serialize_forest(forest);
  if (!predictions.empty()) {
    add_predictions(result, predictions);
  }
  return result;
}

Forest RcppUtilities::deserialize_forest(const Rcpp::List& forest_object) {
  size_t ci_group_size = forest_object["_ci_group_size"];
  size_t num_variables = forest_object["_num_variables"];

  size_t num_trees = forest_object["_num_trees"];
  std::vector<std::unique_ptr<Tree>> trees;
  trees.reserve(num_trees);

  Rcpp::List root_nodes = forest_object["_root_nodes"];
  Rcpp::List child_nodes = forest_object["_child_nodes"];
  Rcpp::List leaf_samples = forest_object["_leaf_samples"];
  Rcpp::List split_vars = forest_object["_split_vars"];
  Rcpp::List split_values = forest_object["_split_values"];
  Rcpp::List drawn_samples = forest_object["_drawn_samples"];
  Rcpp::List send_missing_left = forest_object["_send_missing_left"];

  Rcpp::List prediction_values = forest_object["_pv_values"];
  size_t num_types = forest_object["_pv_num_types"];

  for (size_t t = 0; t < num_trees; t++) {
    trees.emplace_back(new Tree(
                         root_nodes.at(t),
                         child_nodes.at(t),
                         leaf_samples.at(t),
                         split_vars.at(t),
                         split_values.at(t),
                         drawn_samples.at(t),
                         send_missing_left.at(t),
                         PredictionValues(prediction_values.at(t), num_types)));
  }

  return Forest(trees, num_variables, ci_group_size);
}

Rcpp::List RcppUtilities::serialize_forest(Forest& forest) {
  Rcpp::List result;

  result.push_back(forest.get_ci_group_size(), "_ci_group_size");
  result.push_back(forest.get_num_variables(), "_num_variables");

  size_t num_trees = forest.get_trees().size();
  result.push_back(num_trees, "_num_trees");

  Rcpp::List root_nodes(num_trees);
  Rcpp::List child_nodes(num_trees);
  Rcpp::List leaf_samples(num_trees);
  Rcpp::List split_vars(num_trees);
  Rcpp::List split_values(num_trees);
  Rcpp::List drawn_samples(num_trees);
  Rcpp::List send_missing_left(num_trees);
  Rcpp::List prediction_values(num_trees);
  size_t num_types = 0;

  for (size_t t = 0; t < num_trees; t++) {
    // Destructively iterate over the forest by moving the unique_ptr to each tree.
    std::unique_ptr<Tree> tree = std::move(forest.get_trees_().at(t));
    root_nodes[t] = tree->get_root_node();
    child_nodes[t] = tree->get_child_nodes();
    leaf_samples[t] = tree->get_leaf_samples();
    split_vars[t] = tree->get_split_vars();
    split_values[t] = tree->get_split_values();
    drawn_samples[t] = tree->get_drawn_samples();
    send_missing_left[t] = tree->get_send_missing_left();

    prediction_values[t] = tree->get_prediction_values().get_all_values();
    num_types = tree->get_prediction_values().get_num_types();
  }

  result.push_back(root_nodes, "_root_nodes");
  result.push_back(child_nodes, "_child_nodes");
  result.push_back(leaf_samples, "_leaf_samples");
  result.push_back(split_vars, "_split_vars");
  result.push_back(split_values, "_split_values");
  result.push_back(drawn_samples, "_drawn_samples");
  result.push_back(send_missing_left, "_send_missing_left");
  result.push_back(prediction_values, "_pv_values");
  result.push_back(num_types, "_pv_num_types");
  return result;
};

Data RcppUtilities::convert_data(const Rcpp::NumericMatrix& input_data) {
  return Data(input_data.begin(), input_data.nrow(), input_data.ncol());
}

Rcpp::List RcppUtilities::create_prediction_object(const std::vector<Prediction>& predictions) {
  Rcpp::List result;
  add_predictions(result, predictions);
  return result;
};

void RcppUtilities::add_predictions(Rcpp::List& output,
                                    const std::vector<Prediction>& predictions) {
  output.push_back(RcppUtilities::create_prediction_matrix(predictions), "predictions");
  output.push_back(RcppUtilities::create_variance_matrix(predictions), "variance.estimates");
  output.push_back(RcppUtilities::create_error_matrix(predictions), "debiased.error");
  output.push_back(RcppUtilities::create_excess_error_matrix(predictions), "excess.error");
}

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

  Rcpp::NumericMatrix result(predictions.size(), 1);

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

  Rcpp::NumericMatrix result(predictions.size(), 1);

  for (size_t i = 0; i < predictions.size(); i++) {
    const std::vector<double>& error_estimate = predictions[i].get_excess_error_estimates();
    for (size_t j = 0; j < error_estimate.size(); j++) {
      double value = error_estimate[j];
      result(i, j) = value;
    }
  }
  return result;
}
