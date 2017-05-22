/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_FOREST_H_
#define GRF_FOREST_H_

#include <memory>

#include "tree/TreeTrainer.h"
#include "commons/globals.h"
#include "tree/Tree.h"
#include "commons/Data.h"
#include "commons/Observations.h"

class Forest {
public:
  static Forest create(std::vector<std::shared_ptr<Tree>> trees,
                       Data* data,
                       std::unordered_map<size_t, size_t> observables);
  Forest(const std::vector<std::shared_ptr<Tree>>& trees,
         const Observations& observations);

  const Observations& get_observations() const {
    return observations;
  };

  const std::vector<std::shared_ptr<Tree>>& get_trees() const {
    return trees;
  }

protected:
  std::vector<std::shared_ptr<Tree>> trees;
  Observations observations;
};

#endif /* GRF_FOREST_H_ */
