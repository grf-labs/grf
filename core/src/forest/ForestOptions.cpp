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

#include "forest/ForestOptions.h"

ForestOptions::ForestOptions(uint num_trees, uint num_threads, uint random_seed):
    num_trees(num_trees), num_threads(num_threads), random_seed(random_seed) {}

uint ForestOptions::get_num_trees() {
  return num_trees;
}

uint ForestOptions::get_num_threads() {
  return num_threads;
}

uint ForestOptions::get_random_seed() {
  return random_seed;
}
