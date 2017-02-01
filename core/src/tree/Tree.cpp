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

  Authorship: Julie Tibshirani (jtibs@cs.stanford.edu), loosely based on code
  by Marvin Wright (wright@imbs.uni-luebeck.de)
 #-------------------------------------------------------------------------------*/

#include <iterator>
#include "BootstrapSampler.h"

#include "Tree.h"
#include "utility.h"

Tree::Tree(const std::vector<std::vector<size_t>> &child_nodeIDs,
           const std::vector<std::vector<size_t>>& leaf_nodeIDs,
           const std::vector<size_t>& split_varIDs,
           const std::vector<double>& split_values,
           const std::vector<size_t>& oob_sampleIDs,
           const PredictionValues& prediction_values) :
    child_nodeIDs(child_nodeIDs),
    leaf_nodeIDs(leaf_nodeIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    oob_sampleIDs(oob_sampleIDs),
    prediction_values(prediction_values) {}

Tree::~Tree() {}


