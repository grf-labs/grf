/*-------------------------------------------------------------------------------
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

#include <iterator>
#include "sampling/RandomSampler.h"

#include "tree/Tree.h"
#include "commons/utility.h"

namespace grf {

Tree::Tree(size_t root_node,
           const std::vector<std::vector<size_t>>& child_nodes,
           const std::vector<std::vector<size_t>>& leaf_samples,
           const std::vector<size_t>& split_vars,
           const std::vector<double>& split_values,
           const std::vector<size_t>& drawn_samples,
           const PredictionValues& prediction_values) :
    root_node(root_node),
    child_nodes(child_nodes),
    leaf_samples(leaf_samples),
    split_vars(split_vars),
    split_values(split_values),
    drawn_samples(drawn_samples),
    prediction_values(prediction_values) {}

size_t Tree::get_root_node() const {
  return root_node;
}

const std::vector<std::vector<size_t>>& Tree::get_child_nodes() const {
  return child_nodes;
}

const std::vector<std::vector<size_t>>& Tree::get_leaf_samples() const {
  return leaf_samples;
}

const std::vector<size_t>& Tree::get_split_vars() const  {
  return split_vars;
}

const std::vector<double>& Tree::get_split_values() const  {
  return split_values;
}

const std::vector<size_t>& Tree::get_drawn_samples() const  {
  return drawn_samples;
}

const PredictionValues& Tree::get_prediction_values() const  {
  return prediction_values;
}

std::vector<size_t> Tree::find_leaf_nodes(const Data& data,
                                          const std::vector<size_t>& samples) const  {
  std::vector<size_t> prediction_leaf_nodes;
  prediction_leaf_nodes.resize(data.get_num_rows());

  for (size_t sample : samples) {
    size_t node = find_leaf_node(data, sample);
    prediction_leaf_nodes[sample] = node;
  }
  return prediction_leaf_nodes;
}

std::vector<size_t> Tree::find_leaf_nodes(const Data& data,
                                          const std::vector<bool>& valid_samples) const  {
  size_t num_samples = data.get_num_rows();

  std::vector<size_t> prediction_leaf_nodes;
  prediction_leaf_nodes.resize(num_samples);

  for (size_t sample = 0; sample < num_samples; sample++) {
    if (!valid_samples[sample]) {
      continue;
    }

    size_t node = find_leaf_node(data, sample);
    prediction_leaf_nodes[sample] = node;
  }
  return prediction_leaf_nodes;
}

void Tree::set_leaf_samples(const std::vector<std::vector<size_t>>& leaf_samples) {
  this->leaf_samples = leaf_samples;
}

void Tree::set_prediction_values(const PredictionValues& prediction_values) {
  this->prediction_values = prediction_values;
}


size_t Tree::find_leaf_node(const Data& data,
                            size_t sample) const  {
  size_t node = root_node;
  while (true) {
    // Break if terminal node
    if (is_leaf(node)) {
      break;
    }

    // Move to child
    size_t split_var = get_split_vars()[node];
    double value = data.get(sample, split_var);
    if (value <= get_split_values()[node]) {
      // Move to left child
      node = child_nodes[0][node];
    } else {
      // Move to right child
      node = child_nodes[1][node];
    }
  }
  return node;
};

void Tree::prune_empty_leaves() {
  size_t num_nodes = leaf_samples.size();
  for (size_t n = num_nodes; n > root_node; n--) {
    size_t node = n - 1;
    if (is_leaf(node)) {
      continue;
    }

    size_t& left_child = child_nodes[0][node];
    if (!is_leaf(left_child)) {
      prune_node(left_child);
    }

    size_t& right_child = child_nodes[1][node];
    if (!is_leaf(right_child)) {
      prune_node(right_child);
    }
  }
  prune_node(root_node);
}

void Tree::prune_node(size_t& node) {
  size_t left_child = child_nodes[0][node];
  size_t right_child = child_nodes[1][node];

  // If either child is empty, prune this node.
  if (is_empty_leaf(left_child) || is_empty_leaf(right_child)) {
    // Empty out this node.
    child_nodes[0][node] = 0;
    child_nodes[1][node] = 0;

    // If one of the children is not empty, promote it.
    if (!is_empty_leaf(left_child)) {
      node = left_child;
    } else if (!is_empty_leaf(right_child)) {
      node = right_child;
    }
  }
}

bool Tree::is_leaf(size_t node) const  {
  return child_nodes[0][node] == 0 && child_nodes[1][node] == 0;
}

bool Tree::is_empty_leaf(size_t node) const  {
  return is_leaf(node) && leaf_samples[node].empty();
}

} // namespace grf
