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

void ForestTestUtilities::init_default_trainer(ForestTrainer &trainer) {
  init_trainer(trainer, false, 1);
}

void ForestTestUtilities::init_honest_trainer(ForestTrainer& trainer) {
  init_trainer(trainer, true, 1);
}

void ForestTestUtilities::init_trainer(ForestTrainer& trainer,
                                       bool honesty,
                                       uint ci_group_size) {
  uint mtry = 3;
  uint num_trees = 50;
  uint seed = 42;
  uint num_threads = 4;
  uint min_node_size = 1;
  std::set<size_t> no_split_variables;
  std::string split_select_weights_file = "";
  bool sample_with_replacement = true;
  std::string sample_weights_file = "";
  double sample_fraction = ci_group_size > 1 ? 0.35 : 0.7;

  trainer.init(mtry, num_trees, seed, num_threads,
               min_node_size, no_split_variables, split_select_weights_file,
               sample_with_replacement, sample_weights_file, sample_fraction,
               honesty, ci_group_size);
}
