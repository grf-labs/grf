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

#ifndef GRF_SAMPLEWEIGHTCOMPUTER_H
#define GRF_SAMPLEWEIGHTCOMPUTER_H

#include "Eigen/Sparse"
#include "forest/Forest.h"

#include <vector>

namespace grf {

class SampleWeightComputer {
public:
  Eigen::SparseVector<double> compute_weights(size_t sample,
                                              const Forest& forest,
                                              const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
                                              const std::vector<std::vector<bool>>& valid_trees_by_sample,
                                              const Data& train_data) const;
};

} // namespace grf

#endif //GRF_SAMPLEWEIGHTCOMPUTER_H
