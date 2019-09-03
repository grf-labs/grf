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

#include "SplitFrequencyComputer.h"

namespace grf {

std::vector<std::vector<size_t>> SplitFrequencyComputer::compute(const Forest& forest,
                                                                 size_t max_depth) const {
  size_t num_variables = forest.get_num_variables();
  std::vector<std::vector<size_t>> result(max_depth, std::vector<size_t>(num_variables));

  for (const auto& tree : forest.get_trees()) {
    const std::vector<std::vector<size_t>>& child_nodes = tree->get_child_nodes();

    size_t depth = 0;
    std::vector<size_t> level = {tree->get_root_node()};

    while (level.size() > 0 && depth < max_depth) {
      std::vector<size_t> next_level;

      for (size_t node : level) {
        if (tree->is_leaf(node)) {
          continue;
        }

        size_t variable = tree->get_split_vars().at(node);
        result[depth][variable]++;

        next_level.push_back(child_nodes[0][node]);
        next_level.push_back(child_nodes[1][node]);
      }

      level = next_level;
      depth++;
    }
  }
  return result;
}

} // namespace grf
