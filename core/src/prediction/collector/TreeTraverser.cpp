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

#include "TreeTraverser.h"
#include "commons/utility.h"

#include <future>

TreeTraverser::TreeTraverser(uint num_threads) :
    num_threads(num_threads) {}

std::vector<std::vector<size_t>> TreeTraverser::get_leaf_nodes(
    const Forest& forest,
    Data* data,
    bool oob_prediction) const {
  size_t num_trees = forest.get_trees().size();

  std::vector<std::vector<size_t>> leaf_nodes_by_tree;
  leaf_nodes_by_tree.reserve(num_trees);

  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, num_trees - 1, num_threads);

  std::vector<std::future<
      std::vector<std::vector<size_t>>>> futures;
  futures.reserve(thread_ranges.size());

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_trees_batch = thread_ranges[i + 1] - start_index;
    futures.push_back(std::async(std::launch::async,
                                 &TreeTraverser::get_leaf_node_batch,
                                 this,
                                 start_index,
                                 num_trees_batch,
                                 forest,
                                 data,
                                 oob_prediction));
  }

  for (auto& future : futures) {
    std::vector<std::vector<size_t>> leaf_nodes = future.get();
    leaf_nodes_by_tree.insert(leaf_nodes_by_tree.end(),
                              leaf_nodes.begin(),
                              leaf_nodes.end());
  }

  return leaf_nodes_by_tree;
};

std::vector<std::vector<bool>> TreeTraverser::get_valid_trees_by_sample(const Forest& forest,
                                                                        Data* data,
                                                                        bool oob_prediction) const {
  size_t num_trees = forest.get_trees().size();
  size_t num_samples = data->get_num_rows();

  std::vector<std::vector<bool>> result(num_samples, std::vector<bool>(num_trees, true));
  if (oob_prediction) {
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      for (size_t sample : forest.get_trees()[tree_idx]->get_drawn_samples()) {
        result[sample][tree_idx] = false;
      }
    }
  }
  return result;
}

std::vector<std::vector<size_t>> TreeTraverser::get_leaf_node_batch(
    size_t start,
    size_t num_trees,
    const Forest& forest,
    Data* data,
    bool oob_prediction) const {

  size_t num_samples = data->get_num_rows();
  std::vector<std::vector<size_t>> all_leaf_nodes(num_trees);

  for (size_t i = 0; i < num_trees; ++i) {
    std::shared_ptr<Tree> tree = forest.get_trees()[start + i];

    std::vector<bool> valid_samples = get_valid_samples(num_samples, tree, oob_prediction);
    std::vector<size_t> leaf_nodes = tree->find_leaf_nodes(data, valid_samples);
    all_leaf_nodes[i] = leaf_nodes;
  }

  return all_leaf_nodes;
}

std::vector<bool> TreeTraverser::get_valid_samples(size_t num_samples,
                                                     std::shared_ptr<Tree> tree,
                                                     bool oob_prediction) const {
  std::vector<bool> valid_samples(num_samples, true);
  if (oob_prediction) {
    for (size_t sample : tree->get_drawn_samples()) {
      valid_samples[sample] = false;
    }
  }
  return valid_samples;
}
