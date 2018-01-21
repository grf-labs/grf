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

#ifndef GRF_FORESTOPTIONS_H
#define GRF_FORESTOPTIONS_H


#include <tree/TreeOptions.h>
#include <sampling/SamplingOptions.h>
#include "commons/globals.h"

class ForestOptions {
public:
  ForestOptions(uint num_trees,
                uint ci_group_size,
                double sample_fraction,
                uint mtry,
                uint min_node_size,
                bool honesty,
                bool sample_with_replacement,
                uint num_threads,
                uint random_seed);

  uint get_num_trees() const;
  uint get_ci_group_size() const;
  double get_sample_fraction() const;

  const TreeOptions& get_tree_options() const;
  const SamplingOptions& get_sampling_options() const;

  uint get_num_threads() const;
  uint get_random_seed() const;


  uint get_min_node_size() const;

  // TODO(jtibs): check the C++ best practices on mutability, and perhaps
  // replace this with an immutable factory method 'with_min_node_size'.
  void set_min_node_size(uint min_node_size);

private:
  uint num_trees;
  uint ci_group_size;
  double sample_fraction;

  TreeOptions tree_options;
  SamplingOptions sampling_options;

  uint num_threads;
  uint random_seed;
};


#endif //GRF_FORESTOPTIONS_H
