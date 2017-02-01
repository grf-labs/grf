/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_FORESTTESTUTILITIES_H
#define GRADIENTFOREST_FORESTTESTUTILITIES_H

#include "ForestTrainer.h"


class ForestTestUtilities {
public:
  static void init_trainer(ForestTrainer& trainer);
  static void init_honest_trainer(ForestTrainer& trainer);
private:
  static void init_trainer(ForestTrainer& trainer, bool honesty);
};

#endif //GRADIENTFOREST_FORESTTESTUTILITIES_H
