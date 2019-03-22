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

#include "SampleWeightComputer.h"

#include "tree/Tree.h"

std::unordered_map<size_t, double> SampleWeightComputer::compute_weights(size_t sample,
                                                                         const Forest& forest,
                                                                         const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
                                                                         const std::vector<std::vector<bool>>& valid_trees_by_sample) const {
  std::unordered_map<size_t, double> weights_by_sample;

  // Create a list of weighted neighbors for this sample.
  for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
    if (!valid_trees_by_sample[sample][tree_index]) {
      continue;
    }

    const std::vector<size_t>& leaf_nodes = leaf_nodes_by_tree.at(tree_index);
    size_t node = leaf_nodes.at(sample);

    std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
    const std::vector<size_t>& samples = tree->get_leaf_samples()[node];
    if (!samples.empty()) {
      add_sample_weights(samples, weights_by_sample);
    }
  }

  normalize_sample_weights(weights_by_sample);
  return weights_by_sample;
}

void SampleWeightComputer::add_sample_weights(const std::vector<size_t>& samples,
                                              std::unordered_map<size_t, double>& weights_by_sample) const {
  double sample_weight = 1.0 / samples.size();

  for (auto& sample : samples) {
    weights_by_sample[sample] += sample_weight;
  }
}

void SampleWeightComputer::normalize_sample_weights(std::unordered_map<size_t, double>& weights_by_sample) const {
  double total_weight = 0.0;
  for (auto it = weights_by_sample.begin(); it != weights_by_sample.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sample.begin(); it != weights_by_sample.end(); ++it) {
    it->second /= total_weight;
  }
}
