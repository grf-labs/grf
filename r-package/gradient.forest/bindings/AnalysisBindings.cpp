/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <Rcpp.h>
#include <queue>
#include <vector>

#include "analysis/SplitFrequencyComputer.h"
#include "commons/globals.h"
#include "forest/Forest.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_split_frequencies(Rcpp::List forest_object,
                                              size_t max_depth) {
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  SplitFrequencyComputer computer;
  std::vector<std::vector<size_t>> variable_frequencies = computer.compute(forest, max_depth);

  size_t num_variables = forest.get_num_variables();
  Rcpp::NumericMatrix result(num_variables, max_depth);
  for (size_t var = 0; var < num_variables; var++) {
    const std::vector<size_t>& frequencies = variable_frequencies.at(var);
    for (size_t depth = 0; depth < frequencies.size(); depth++) {
      double frequency = frequencies[depth];
      result(var, depth) = frequency;
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List examine_tree(Rcpp::List forest_object,
                        size_t tree_index) {
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

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
  size_t node_index = 2; // We start at 2 because R is 1-indexed.

  std::vector<Rcpp::List> node_objects;
  while (frontier.size() > 0) {
    size_t node = frontier.front();
    Rcpp::List node_object;

    if (tree->is_leaf(node)) {
      node_object.push_back(true, "is_leaf");
      node_object.push_back(leaf_samples.at(node), "samples");
    } else {
      node_object.push_back(false, "is_leaf");
      node_object.push_back(split_vars.at(node), "split_variable");
      node_object.push_back(split_values.at(node), "split_value");

      node_object.push_back(node_index++, "left_child");
      frontier.push(child_nodes[0][node]);

      node_object.push_back(node_index++, "right_child");
      frontier.push(child_nodes[1][node]);
    }

    frontier.pop();
    node_objects.push_back(node_object);
  }

  Rcpp::List result;
  result.push_back(tree->get_oob_samples(), "oob_samples");
  result.push_back(node_objects, "nodes");
  return result;
}
