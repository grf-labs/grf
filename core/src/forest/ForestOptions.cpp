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
 
#include <thread>
#include <stdexcept>

#include "forest/ForestOptions.h"
#include "tree/TreeOptions.h"

namespace grf {

ForestOptions::ForestOptions(uint num_trees,
                             size_t ci_group_size,
                             double sample_fraction,
                             uint mtry,
                             uint min_node_size,
                             bool honesty,
                             double honesty_fraction,
                             bool honesty_prune_leaves,
                             double alpha,
                             double imbalance_penalty,
                             uint num_threads,
                             uint random_seed,
                             const std::vector<size_t>& sample_clusters,
                             uint samples_per_cluster):
    ci_group_size(ci_group_size),
    sample_fraction(sample_fraction),
    tree_options(mtry, min_node_size, honesty, honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty),
    sampling_options(samples_per_cluster, sample_clusters),
    random_seed(random_seed) {

  this->num_threads = validate_num_threads(num_threads);

  // If necessary, round the number of trees up to a multiple of
  // the confidence interval group size.
  this->num_trees = num_trees + (num_trees % ci_group_size);

  if (ci_group_size > 1 && sample_fraction > 0.5) {
    throw std::runtime_error("When confidence intervals are enabled, the"
        " sampling fraction must be less than 0.5.");
  }
}

uint ForestOptions::get_num_trees() const {
  return num_trees;
}

size_t ForestOptions::get_ci_group_size() const {
  return ci_group_size;
}

double ForestOptions::get_sample_fraction() const {
  return sample_fraction;
}

const TreeOptions& ForestOptions::get_tree_options() const {
  return tree_options;
}

const SamplingOptions& ForestOptions::get_sampling_options() const {
  return sampling_options;
}

uint ForestOptions::get_num_threads() const {
  return num_threads;
}

uint ForestOptions::get_random_seed() const {
  return random_seed;
}

uint ForestOptions::validate_num_threads(uint num_threads) {
  if (num_threads == DEFAULT_NUM_THREADS) {
    return std::thread::hardware_concurrency();
  } else if (num_threads > 0) {
    return num_threads;
  } else {
    throw std::runtime_error("A negative number of threads was provided.");
  }
}

} // namespace grf
