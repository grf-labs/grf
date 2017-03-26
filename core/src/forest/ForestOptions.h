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

#ifndef GRADIENTFOREST_FORESTOPTIONS_H
#define GRADIENTFOREST_FORESTOPTIONS_H


#include "commons/globals.h"

class ForestOptions {
public:
  ForestOptions(uint num_trees, uint num_threads, uint random_seed);

  uint get_num_trees();
  uint get_num_threads();
  uint get_random_seed();

private:
  uint num_trees;
  uint num_threads;
  uint random_seed;
};


#endif //GRADIENTFOREST_FORESTOPTIONS_H
