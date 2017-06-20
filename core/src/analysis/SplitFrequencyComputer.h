/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#ifndef GRF_VARIABLEIMPORTANCECOMPUTER_H
#define GRF_VARIABLEIMPORTANCECOMPUTER_H


#include "forest/Forest.h"

/**
 * Computes a matrix of split depth by variable ID, where each value is
 * the number of times the variable was split on at that depth.
 *
 * forest: the forest for which split frequencies should be computed
 * max_depth: the maximum depth of splits to consider, exclusive
 */
class SplitFrequencyComputer {
public:
  std::vector<std::vector<size_t>> compute(const Forest& forest,
                                           size_t max_depth);
};


#endif //GRF_VARIABLEIMPORTANCECOMPUTER_H
