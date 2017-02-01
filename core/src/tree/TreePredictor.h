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

  Author: Julie Tibshirani (jtibs@cs.stanford.edu)
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_TREEPREDICTOR_H
#define GRADIENTFOREST_TREEPREDICTOR_H

#include <memory>

#include "Data.h"
#include "Tree.h"

class TreePredictor {
public:
  std::vector<size_t> get_terminal_nodeIDs(std::shared_ptr<Tree>,
                                           Data* prediction_data,
                                           const std::vector<size_t>& sampleIDs);
};


#endif //GRADIENTFOREST_TREEPREDICTOR_H
