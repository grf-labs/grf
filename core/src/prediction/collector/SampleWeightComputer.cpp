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

#include <algorithm>
#include <array>
#include "SampleWeightComputer.h"
#include "tree/Tree.h"

std::unordered_map<size_t, double> SampleWeightComputer::compute_weights(size_t sample,
                                                                         const Forest& forest,
                                                                         const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
                                                                         const std::vector<std::vector<bool>>& valid_trees_by_sample) const {
  size_t n = forest.get_observations().get_num_samples();

  // Temporary vector of unnormalized neighbor weights. Later, zero-weighted neighbors will be discarded.
  std::vector<double> raw_weights(n, 0.);
  double total_weights = 0.;
  
  // Loops through each leaf node containing 'sample', and stores who are the neighbors of 'sample'.
  for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
    if (!valid_trees_by_sample[sample][tree_index]) {
      continue;
    }

    const std::vector<size_t>& leaf_nodes = leaf_nodes_by_tree.at(tree_index);
    size_t node = leaf_nodes.at(sample);

    std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
    const std::vector<size_t>& samples = tree->get_leaf_samples()[node];
    
    for (auto&s : samples) {
      ++raw_weights[s];
      ++total_weights;
    }
  }
  
  // This map will now hold only neighbors with non-zero weights.
  std::unordered_map<size_t, double> weights_by_sample;
  for (size_t i=0; i < n; ++i) {
      double w = raw_weights[i];
      if (w > 0) {
        weights_by_sample[i] = w / total_weights;
      }
  }
  
  return weights_by_sample;
}
