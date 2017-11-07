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

#ifndef GRF_FORESTTESTUTILITIES_H
#define GRF_FORESTTESTUTILITIES_H

#include "forest/ForestTrainer.h"


class ForestTestUtilities {
public:
    static void init_default_trainer(ForestTrainer &trainer);
    static void init_honest_trainer(ForestTrainer& trainer);
    static void init_trainer(ForestTrainer& trainer,
                             bool honesty,
                             uint min_node_size);
    static void init_trainer_test(ForestTrainer& trainer,
                           bool honesty,
                           uint ci_group_size,
                           uint number_of_tree,
                           const double *query,
                           uint variable_dimensions);
};

#endif //GRF_FORESTTESTUTILITIES_H
