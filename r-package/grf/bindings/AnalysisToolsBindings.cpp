/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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
#include <queue>
#include <vector>

#include "analysis/SplitFrequencyComputer.h"
#include "commons/globals.h"
#include "forest/Forest.h"
#include "prediction/collector/SampleWeightComputer.h"
#include "prediction/collector/TreeTraverser.h"

#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_split_frequencies(Rcpp::List forest_object,
                                              size_t max_depth) {
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  SplitFrequencyComputer computer;
  std::vector<std::vector<size_t>> split_frequencies = computer.compute(forest, max_depth);

  size_t num_variables = forest.get_num_variables();
  Rcpp::NumericMatrix result(max_depth, num_variables);
  for (size_t depth = 0; depth < split_frequencies.size(); depth++) {
    const std::vector<size_t>& frequencies = split_frequencies.at(depth);
    for (size_t var = 0; var < num_variables; var++) {
      double frequency = frequencies[var];
        result(depth, var) = frequency;
      }
    }
  return result;
}

Eigen::SparseMatrix<double> compute_sample_weights(Rcpp::List forest_object,
                                                   Rcpp::NumericMatrix input_data,
                                                   Eigen::SparseMatrix<double> sparse_input_data,
                                                   unsigned int num_threads,
                                                   bool oob_prediction) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);
  num_threads = ForestOptions::validate_num_threads(num_threads);

  TreeTraverser tree_traverser(num_threads);
  SampleWeightComputer weight_computer;

  std::vector<std::vector<size_t>> leaf_nodes_by_tree = tree_traverser.get_leaf_nodes(forest, data, oob_prediction);
  std::vector<std::vector<bool>> trees_by_sample = tree_traverser.get_valid_trees_by_sample(forest, data, oob_prediction);

  size_t num_samples = data->get_num_rows();
  size_t num_neighbors = forest.get_observations().get_num_samples();

  Eigen::SparseMatrix<double> result(num_samples, num_neighbors);

  for (size_t sample = 0; sample < num_samples; sample++) {
    std::unordered_map<size_t, double> weights = weight_computer.compute_weights(
        sample, forest, leaf_nodes_by_tree, trees_by_sample);
    for (auto it = weights.begin(); it != weights.end(); it++) {
      size_t neighbor = it->first;
      double weight = it->second;
      result.insert(sample, neighbor) = weight;
    }
  }
  result.makeCompressed();
  return result;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_weights(Rcpp::List forest_object,
                                            Rcpp::NumericMatrix input_data,
                                            Eigen::SparseMatrix<double> sparse_input_data,
                                            unsigned int num_threads) {
  return compute_sample_weights(forest_object, input_data, sparse_input_data, num_threads, false);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_weights_oob(Rcpp::List forest_object,
                                                Rcpp::NumericMatrix input_data,
                                                Eigen::SparseMatrix<double> sparse_input_data,
                                                unsigned int num_threads) {
  return compute_sample_weights(forest_object, input_data, sparse_input_data, num_threads, true);
}

// [[Rcpp::export]]
Rcpp::List deserialize_tree(Rcpp::List forest_object,
                            size_t tree_index) {
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  tree_index--; // Decrement since R is one-indexed.
  size_t num_trees = forest.get_trees().size();
  if (tree_index >= num_trees) {
    throw std::runtime_error("The provided tree index is not valid.");
  }

  std::shared_ptr<Tree> tree = forest.get_trees().at(tree_index);
  const std::vector<std::vector<size_t>>& child_nodes = tree->get_child_nodes();
  const std::vector<std::vector<size_t>>& leaf_samples = tree->get_leaf_samples();

  const std::vector<size_t>& split_vars = tree->get_split_vars();
  const std::vector<double>& split_values = tree->get_split_values();

  std::queue<size_t> frontier;
  frontier.push(tree->get_root_node());
  size_t node_index = 1;

  std::vector<Rcpp::List> node_objects;

  // Note that since R is 1-indexed, we add '1' below to array indices.
  while (frontier.size() > 0) {
    size_t node = frontier.front();
    Rcpp::List node_object;

    if (tree->is_leaf(node)) {
      node_object.push_back(true, "is_leaf");
      node_object.push_back(leaf_samples.at(node), "samples");
    } else {
      node_object.push_back(false, "is_leaf");
      node_object.push_back(split_vars.at(node) + 1, "split_variable"); // R is 1-indexed.
      node_object.push_back(split_values.at(node), "split_value");

      node_object.push_back(node_index + 1, "left_child");
      frontier.push(child_nodes[0][node]);
      node_index++;

      node_object.push_back(node_index + 1, "right_child");
      frontier.push(child_nodes[1][node]);
      node_index++;
    }

    frontier.pop();
    node_objects.push_back(node_object);
  }

  const std::vector<size_t>& drawn_samples = tree->get_drawn_samples();

  Rcpp::List result;
  result.push_back(drawn_samples.size(), "num_samples");
  result.push_back(drawn_samples, "drawn_samples");
  result.push_back(node_objects, "nodes");
  return result;
}
