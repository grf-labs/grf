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

#include "commons/DefaultData.h"
#include "forest/Forest.h"

Forest Forest::create(std::vector<std::shared_ptr<Tree>> trees,
                      const ForestOptions& forest_options,
                      const Data* data) {
  size_t num_independent_variables = data->get_num_cols() -
      data->get_disallowed_split_variables().size();
  return Forest(trees, num_independent_variables, forest_options.get_ci_group_size());
}

Forest::Forest(const std::vector<std::shared_ptr<Tree>>& trees,
               size_t num_variables,
               size_t ci_group_size):
  trees(trees),
  num_variables(num_variables),
  ci_group_size(ci_group_size) {}

Forest Forest::merge(const std::vector<std::shared_ptr<Forest>>& forests) {

  std::vector<std::shared_ptr<Tree>> all_trees;
  const size_t num_variables = forests[0]->get_num_variables();
  const size_t ci_group_size = forests[0]->get_ci_group_size();

  for (auto& forest : forests) {
    auto& trees = forest->get_trees();
    all_trees.insert(all_trees.end(), trees.begin(), trees.end());

    if (forest->get_ci_group_size() != ci_group_size) {
      throw std::runtime_error("All forests being merged must have the same ci_group_size.");
    }
  }

  return Forest(all_trees, num_variables, ci_group_size);
}

const std::vector<std::shared_ptr<Tree>>& Forest::get_trees() const {
  return trees;
}

const size_t Forest::get_num_variables() const {
  return num_variables;
}

const size_t Forest::get_ci_group_size() const {
  return ci_group_size;
}
