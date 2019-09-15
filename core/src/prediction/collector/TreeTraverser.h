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

#ifndef GRF_TREETRAVERSER_H
#define GRF_TREETRAVERSER_H

#include "forest/Forest.h"

namespace grf {

class TreeTraverser {
public:
  TreeTraverser(uint num_threads);

  std::vector<std::vector<size_t>> get_leaf_nodes(
      const Forest& forest,
      const Data& data,
      bool oob_prediction) const;

  std::vector<std::vector<bool>> get_valid_trees_by_sample(const Forest& forest,
                                                           const Data& data,
                                                           bool oob_prediction) const;

private:
  std::vector<std::vector<size_t>> get_leaf_node_batch(
      size_t start,
      size_t num_trees,
      const Forest& forest,
      const Data& data,
      bool oob_prediction) const;

  std::vector<bool> get_valid_samples(size_t num_samples,
                                      const std::unique_ptr<Tree>& tree,
                                      bool oob_prediction) const;

  uint num_threads;
};

} // namespace grf

#endif //GRF_TREETRAVERSER_H
