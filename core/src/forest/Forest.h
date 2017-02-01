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
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_FOREST_H_
#define GRADIENTFOREST_FOREST_H_

#include <memory>

#include "TreeTrainer.h"
#include "globals.h"
#include "Tree.h"
#include "Data.h"
#include "PredictionStrategy.h"
#include "Observations.h"

class Forest {
public:
  static Forest create(std::vector<std::shared_ptr<Tree>> trees,
                       Data* data,
                       std::unordered_map<std::string, size_t> observables);
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

#endif /* GRADIENTFOREST_FOREST_H_ */
