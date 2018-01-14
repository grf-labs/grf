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

#include "utilities/ForestTestUtilities.h"
#include "forest/ForestTrainer.h"

ForestOptions ForestTestUtilities::default_options() {
  return default_options(false, 1);
}

ForestOptions ForestTestUtilities::default_honest_options() {
  return default_options(true, 1);
}

ForestOptions ForestTestUtilities::default_options(bool honesty,
                                                   uint ci_group_size) {
  uint mtry = 3;
  uint num_trees = 50;
  uint seed = 42;
  uint num_threads = 4;
  uint min_node_size = 1;
  bool sample_with_replacement = true;
  double sample_fraction = ci_group_size > 1 ? 0.35 : 0.7;

  return ForestOptions(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
                       honesty, sample_with_replacement, num_threads, seed);
}
