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

class Forest {
public:
  static Forest create(std::vector<std::shared_ptr<Tree>> trees,
                       const ForestOptions& forest_options,
                       const Data* data);

  Forest(const std::vector<std::shared_ptr<Tree>>& trees,
         size_t num_variables,
         size_t ci_group_size);

  const std::vector<std::shared_ptr<Tree>>& get_trees() const;

  const size_t get_num_variables() const;
  const size_t get_ci_group_size() const;

  static Forest merge(const std::vector<std::shared_ptr<Forest>>& forests);
  
private:
  std::vector<std::shared_ptr<Tree>> trees;
  size_t num_variables;
  size_t ci_group_size;
};

#endif /* GRF_FOREST_H_ */
