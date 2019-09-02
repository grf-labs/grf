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

#ifndef GRF_TREE_H_
#define GRF_TREE_H_

#include <vector>
#include <iostream>

#include "commons/globals.h"
#include "commons/DefaultData.h"
#include "sampling/RandomSampler.h"
#include "prediction/PredictionValues.h"
#include "splitting/SplittingRule.h"

namespace grf {

class Tree {
public:
  Tree(size_t root_node,
       const std::vector<std::vector<size_t>>& child_nodes,
       const std::vector<std::vector<size_t>>& leaf_samples,
       const std::vector<size_t>& split_vars,
       const std::vector<double>& split_values,
       const std::vector<size_t>& drawn_samples,
       const PredictionValues& prediction_values);

  /**
   * Given test data and a list of sample IDs, recurses down the tree to find
   * the leaf node IDs that those samples belong in.
   *
   * @param data: the data matrix containing all test samples.
   * @param samples: a list of sample IDs whose leaf nodes should be calculated.
   * @return The resulting node IDs for each sample ID in the input. For efficiency reasons, this
   * vector's length will be equal to the total number of test samples, and the index for each
   * requested sample ID will contain the corresponding node ID. All other values will be 0.
   */
  std::vector<size_t> find_leaf_nodes(const Data& data,
                                      const std::vector<size_t>& samples) const;

  /**
   * Given test data and a vector indicating which samples to consider, recurses
   * down the tree to find the leaf node IDs that those samples belong in.
   *
   * @param data: the data matrix containing all test samples.
   * @param samples: a vector indicating which samples should have their leaf nodes
   * calculated. This vector's length should be equal to the total number of samples, and
   * contain 'true' for a sample ID's index if that sample's leaf node should be calculated.
   * @return The resulting node IDs for each sample ID in the input. For efficiency reasons, this
   * vector's length will be equal to the total number of test samples, and the index for each
   * requested sample ID will contain the corresponding node ID. All other values will be 0.
   */
  std::vector<size_t> find_leaf_nodes(const Data& data,
                                      const std::vector<bool>& valid_samples) const;
  /**
   * Removes all empty leaf nodes.
   *
   * When re-populating the leaves of an honest tree, certain leaf nodes may become empty.
   * This procedure prunes those nodes, so that each node is either a non-empty leaf, or
   * has two non-empty subtrees for children.
   */
  void prune_empty_leaves();

  /**
   * The ID of the root node for this tree. Note that this is usually 0, but may not always
   * be as the top of the tree can be pruned.
   */
  size_t get_root_node() const;

  /**
   * A vector containing two vectors: the first gives the ID of the left child for every
   * node, and the second gives the ID of the right child. If a node is a leaf, the entries
   * for both the left and right children will be '0'.
   */
  const std::vector<std::vector<size_t>>& get_child_nodes() const;

  /**
   * Specifies the samples that each node contains. Note that only leaf nodes will contain
   * a non-empty vector of sample IDs.
   */
  const std::vector<std::vector<size_t>>& get_leaf_samples() const;

  /**
   * For each split, the ID of the variable that was chosen to split on.
   */
  const std::vector<size_t>& get_split_vars() const;

  /**
   * For each split, the value of the variable that was chosen to split on.
   */
  const std::vector<double>& get_split_values() const;

  /**
   * The sample IDs that were not drawn in creating this tree. For honest trees,
   * this excludes both samples that went into growing the tree, as well as samples
   * used to repopulate the leaves.
   */
  const std::vector<size_t>& get_drawn_samples() const;

  /**
   * Optional summary values about the samples in each leaf. Note that this will only
   * be non-empty if the tree was trained with an 'optimized' prediction strategy.
   */
  const PredictionValues& get_prediction_values() const;

  /**
   * Given a node ID, returns true if the node represents a leaf in this tree (in
   * particular, the node has no children).
   */
  bool is_leaf(size_t node) const;

  /**
   * Sets the contents of this tree's leaf nodes. Please see
   * Tree::get_leaf_samples for a description of this variable.
   */
  void set_leaf_samples(const std::vector<std::vector<size_t>>& leaf_samples);

  /**
   * Sets the contents of this tree's prediction values. Please see
   * Tree::get_prediction_values for a description of this variable.
   */
  void set_prediction_values(const PredictionValues& prediction_values);

private:
  size_t find_leaf_node(const Data& data,
                        size_t sample) const;
  void prune_node(size_t& node);
  bool is_empty_leaf(size_t node) const;

  size_t root_node;
  std::vector<std::vector<size_t>> child_nodes;
  std::vector<std::vector<size_t>> leaf_samples;
  std::vector<size_t> split_vars;
  std::vector<double> split_values;

  std::vector<size_t> drawn_samples;

  PredictionValues prediction_values;
};

} // namespace grf

#endif /* GRF_TREE_H_ */
