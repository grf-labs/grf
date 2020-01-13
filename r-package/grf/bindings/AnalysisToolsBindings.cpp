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

using namespace grf;

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_split_frequencies(SEXP forest_xptr,
                                              size_t max_depth) {
  Rcpp::XPtr<Forest> forest(forest_xptr);

  SplitFrequencyComputer computer;
  std::vector<std::vector<size_t>> split_frequencies = computer.compute(*forest, max_depth);

  size_t num_variables = forest->get_num_variables();
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

Eigen::SparseMatrix<double> compute_sample_weights(SEXP forest_xptr,
                                                   Rcpp::NumericMatrix train_matrix,
                                                   Eigen::SparseMatrix<double> sparse_train_matrix,
                                                   Rcpp::NumericMatrix test_matrix,
                                                   Eigen::SparseMatrix<double> sparse_test_matrix,
                                                   unsigned int num_threads,
                                                   bool oob_prediction) {
  std::unique_ptr<Data> train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  std::unique_ptr<Data> data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);
  Rcpp::XPtr<Forest> forest(forest_xptr);
  num_threads = ForestOptions::validate_num_threads(num_threads);

  TreeTraverser tree_traverser(num_threads);
  SampleWeightComputer weight_computer;

  std::vector<std::vector<size_t>> leaf_nodes_by_tree = tree_traverser.get_leaf_nodes(*forest, *data, oob_prediction);
  std::vector<std::vector<bool>> trees_by_sample = tree_traverser.get_valid_trees_by_sample(*forest, *data, oob_prediction);

  size_t num_samples = data->get_num_rows();
  size_t num_neighbors = train_data->get_num_rows();

  // From http://eigen.tuxfamily.org/dox/group__TutorialSparse.html:
  // Filling a sparse matrix effectively
  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(num_neighbors);
  Eigen::SparseMatrix<double> result(num_samples, num_neighbors);

  for (size_t sample = 0; sample < num_samples; sample++) {
    std::unordered_map<size_t, double> weights = weight_computer.compute_weights(
        sample, *forest, leaf_nodes_by_tree, trees_by_sample);
    for (auto it = weights.begin(); it != weights.end(); it++) {
      size_t neighbor = it->first;
      double weight = it->second;
      triplet_list.emplace_back(sample, neighbor, weight);
    }
  }
  result.setFromTriplets(triplet_list.begin(), triplet_list.end());

  return result;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_weights(SEXP forest_xptr,
                                            Rcpp::NumericMatrix train_matrix,
                                            Eigen::SparseMatrix<double> sparse_train_matrix,
                                            Rcpp::NumericMatrix test_matrix,
                                            Eigen::SparseMatrix<double> sparse_test_matrix,
                                            unsigned int num_threads) {
  return compute_sample_weights(forest_xptr, train_matrix, sparse_test_matrix,
                                test_matrix, sparse_test_matrix, num_threads, false);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_weights_oob(SEXP forest_xptr,
                                                Rcpp::NumericMatrix test_matrix,
                                                Eigen::SparseMatrix<double> sparse_test_matrix,
                                                unsigned int num_threads) {
  return compute_sample_weights(forest_xptr, test_matrix, sparse_test_matrix,
                                test_matrix, sparse_test_matrix, num_threads, true);
}

// [[Rcpp::export]]
Rcpp::List deserialize_tree(SEXP forest_xptr,
                            size_t tree_index) {
  Rcpp::XPtr<Forest> forest(forest_xptr);

  tree_index--; // Decrement since R is one-indexed.
  size_t num_trees = forest->get_trees().size();
  if (tree_index < 0 || tree_index >= num_trees) {
    throw std::runtime_error("The provided tree index is not valid.");
  }

  const std::unique_ptr<Tree>& tree = forest->get_trees().at(tree_index);
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

      std::vector<size_t> samples;
      samples.reserve(leaf_samples.at(node).size());
      for (size_t index : leaf_samples.at(node)) {
        samples.push_back(index + 1); // R is 1-indexed.
      }
      node_object.push_back(samples, "samples");
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

  std::vector<size_t> drawn_samples(tree->get_drawn_samples());
  for (size_t& index : drawn_samples) index += 1; //R is 1-indexed.


  Rcpp::List result;
  result.push_back(drawn_samples.size(), "num_samples");
  result.push_back(drawn_samples, "drawn_samples");
  result.push_back(node_objects, "nodes");
  return result;
}

// [[Rcpp::export]]
Rcpp::List merge(const Rcpp::List forest_xptr_list) {
  std::vector<std::unique_ptr<Tree>> all_trees;

  auto forest = Rcpp::as< Rcpp::XPtr<Forest> > (forest_xptr_list.at(0));
  const size_t num_variables = forest->get_num_variables();
  const size_t ci_group_size = forest->get_ci_group_size();

  for (auto& opaque_forest : forest_xptr_list) {
    auto forest = Rcpp::as< Rcpp::XPtr<Forest> > (opaque_forest);
    auto& trees = forest->get_trees_();
    all_trees.insert(all_trees.end(),
                     std::make_move_iterator(trees.begin()),
                     std::make_move_iterator(trees.end()));
    if (forest->get_ci_group_size() != ci_group_size) {
      throw std::runtime_error("All forests being merged must have the same ci_group_size.");
    }
  }

  Forest big_forest = Forest(all_trees, num_variables, ci_group_size);
  return RcppUtilities::create_forest_object(big_forest);
}

// [[Rcpp::export]]
Rcpp::List internal_save(SEXP forest_xptr) {
  return RcppUtilities::serialize_forest(forest_xptr);
}

// [[Rcpp::export]]
Rcpp::List internal_load(const Rcpp::List& forest_object) {
  Forest forest = RcppUtilities::deserialize_forest(forest_object);
  return RcppUtilities::create_forest_object(forest);
}
