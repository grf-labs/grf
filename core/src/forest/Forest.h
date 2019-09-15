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

#ifndef GRF_FOREST_H_
#define GRF_FOREST_H_

#include "commons/DefaultData.h"
#include "commons/globals.h"
#include "forest/ForestOptions.h"
#include "tree/TreeTrainer.h"
#include "tree/Tree.h"

namespace grf {

class Forest {
public:
  Forest(std::vector<std::unique_ptr<Tree>>& trees,
         size_t num_variables,
         size_t ci_group_size);

  Forest(Forest&& forest);

  const std::vector<std::unique_ptr<Tree>>& get_trees() const;

  /**
   * A method intended for internal use that allows the list of
   * trees to be modified.
   */
  std::vector<std::unique_ptr<Tree>>& get_trees_();

  const size_t get_num_variables() const;
  const size_t get_ci_group_size() const;

  /**
   * Merges the given forests into a single forest. The new forest
   * will contain all the trees from the smaller forests.
   *
   * NOTE: this is a destructive operation -- the original forests cannot
   * be used after they are merged together.
   */
  static Forest merge(std::vector<Forest>& forests);
  
private:
  std::vector<std::unique_ptr<Tree>> trees;
  size_t num_variables;
  size_t ci_group_size;
  DISALLOW_COPY_AND_ASSIGN(Forest);
};

} // namespace grf

#endif /* GRF_FOREST_H_ */
