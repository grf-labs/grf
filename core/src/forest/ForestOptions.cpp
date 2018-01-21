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
#include <random>
#include "forest/ForestOptions.h"

ForestOptions::ForestOptions(uint num_trees,
                             uint ci_group_size,
                             double sample_fraction,
                             uint mtry,
                             uint min_node_size,
                             bool honesty,
                             bool sample_with_replacement,
                             uint num_threads,
                             uint random_seed):
    ci_group_size(ci_group_size),
    sample_fraction(sample_fraction),
    tree_options(mtry, min_node_size, honesty),
    sampling_options(sample_with_replacement) {

  if (num_threads == DEFAULT_NUM_THREADS) {
    this->num_threads = std::thread::hardware_concurrency();
  } else {
    this->num_threads = num_threads;
  }

  // If necessary, round the number of trees up to a multiple of
  // the confidence interval group size.
  this->num_trees = num_trees + (num_trees % ci_group_size);

  if (ci_group_size > 1 && sample_fraction > 0.5) {
    throw std::runtime_error("When confidence intervals are enabled, the"
        " sampling fraction must be less than 0.5.");
  }

  if (random_seed != 0) {
    this->random_seed = random_seed;
  } else {
    std::random_device random_device;
    this->random_seed = random_device();
  }
}

uint ForestOptions::get_num_trees() const {
  return num_trees;
}

uint ForestOptions::get_ci_group_size() const {
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

uint ForestOptions::get_min_node_size() const {
  return tree_options.get_min_node_size();
}

void ForestOptions::set_min_node_size(uint min_node_size) {
  return tree_options.set_min_node_size(min_node_size);
}
