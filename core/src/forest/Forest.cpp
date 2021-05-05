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

#include <stdexcept>

#include "commons/Data.h"
#include "forest/Forest.h"

namespace grf {

Forest::Forest(std::vector<std::unique_ptr<Tree>>& trees,
               size_t num_variables,
               size_t ci_group_size) {
  this->trees.insert(this->trees.end(),
                     std::make_move_iterator(trees.begin()),
                     std::make_move_iterator(trees.end()));
  this->num_variables = num_variables;
  this->ci_group_size = ci_group_size;
}

Forest::Forest(Forest&& forest) {
  this->trees.insert(this->trees.end(),
                     std::make_move_iterator(forest.trees.begin()),
                     std::make_move_iterator(forest.trees.end()));
  this->num_variables = forest.num_variables;
  this->ci_group_size = forest.ci_group_size;
}

Forest Forest::merge(std::vector<Forest>& forests) {
  std::vector<std::unique_ptr<Tree>> all_trees;
  const size_t num_variables = forests.at(0).get_num_variables();
  const size_t ci_group_size = forests.at(0).get_ci_group_size();

  for (auto& forest : forests) {
    auto& trees = forest.get_trees_();
    all_trees.insert(all_trees.end(),
                     std::make_move_iterator(trees.begin()),
                     std::make_move_iterator(trees.end()));

    if (forest.get_ci_group_size() != ci_group_size) {
      throw std::runtime_error("All forests being merged must have the same ci_group_size.");
    }
  }

  return Forest(all_trees, num_variables, ci_group_size);
}

const std::vector<std::unique_ptr<Tree>>& Forest::get_trees() const {
  return trees;
}

std::vector<std::unique_ptr<Tree>>& Forest::get_trees_() {
  return trees;
}

const size_t Forest::get_num_variables() const {
  return num_variables;
}

const size_t Forest::get_ci_group_size() const {
  return ci_group_size;
}

} // namespace grf
