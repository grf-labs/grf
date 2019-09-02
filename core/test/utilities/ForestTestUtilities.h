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

using namespace grf;

class ForestTestUtilities {
public:
  static ForestOptions default_options();
  static ForestOptions default_honest_options();

  static ForestOptions default_options(bool honesty, size_t ci_group_size);
};

#endif //GRF_FORESTTESTUTILITIES_H
