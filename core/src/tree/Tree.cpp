/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <iterator>
#include "sampling/BootstrapSampler.h"

#include "tree/Tree.h"
#include "commons/utility.h"

Tree::Tree(size_t root_nodeID,
           const std::vector<std::vector<size_t>>& child_nodeIDs,
           const std::vector<std::vector<size_t>>& leaf_nodeIDs,
           const std::vector<size_t>& split_varIDs,
           const std::vector<double>& split_values,
           const std::vector<size_t>& oob_sampleIDs,
           const PredictionValues& prediction_values) :
    root_nodeID(root_nodeID),
    child_nodeIDs(child_nodeIDs),
    leaf_nodeIDs(leaf_nodeIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    oob_sampleIDs(oob_sampleIDs),
    prediction_values(prediction_values) {}

std::vector<size_t> Tree::find_leaf_nodeIDs(Data* prediction_data,
                                            const std::vector<size_t> &sampleIDs) {
  bool use_subsample = !sampleIDs.empty();
  const std::vector<std::vector<size_t>>& child_nodeIDs = get_child_nodeIDs();

  std::vector<size_t> prediction_leaf_nodeIDs;
  prediction_leaf_nodeIDs.resize(prediction_data->get_num_rows());

  size_t num_samples_predict = use_subsample ? sampleIDs.size() : prediction_data->get_num_rows();

  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sampleID = use_subsample ? sampleIDs[i] : i;

    size_t nodeID = root_nodeID;
    while (true) {
      // Break if terminal node
      if (is_leaf(nodeID)) {
        break;
      }

      // Move to child
      size_t split_varID = get_split_varIDs()[nodeID];
      double value = prediction_data->get(sampleID, split_varID);
      if (value <= get_split_values()[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

    prediction_leaf_nodeIDs[sampleID] = nodeID;
  }
  return prediction_leaf_nodeIDs;
}

void Tree::prune_empty_leaves() {
  size_t num_nodes = leaf_nodeIDs.size();
  for (size_t n = num_nodes; n > root_nodeID; n--) {
    size_t node = n - 1;
    if (is_leaf(node)) {
      continue;
    }

    size_t& left_child = child_nodeIDs[0][node];
    if (!is_leaf(left_child)) {
      prune_node(left_child);
    }

    size_t& right_child = child_nodeIDs[1][node];
    if (!is_leaf(right_child)) {
      prune_node(right_child);
    }
  }
  prune_node(root_nodeID);
}

void Tree::prune_node(size_t& node) {
  size_t left_child = child_nodeIDs[0][node];
  size_t right_child = child_nodeIDs[1][node];

  // If either child is empty, prune this node.
  if (is_empty_leaf(left_child) || is_empty_leaf(right_child)) {
    // Empty out this node.
    child_nodeIDs[0][node] = 0;
    child_nodeIDs[1][node] = 0;

    // If one of the children is not empty, promote it.
    if (!is_empty_leaf(left_child)) {
      node = left_child;
    } else if (!is_empty_leaf(right_child)) {
      node = right_child;
    }
  }
}

bool Tree::is_leaf(size_t node) {
  return child_nodeIDs[0][node] == 0 && child_nodeIDs[1][node] == 0;
}

bool Tree::is_empty_leaf(size_t node) {
  return is_leaf(node) && leaf_nodeIDs[node].empty();
}

