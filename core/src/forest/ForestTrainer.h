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

#ifndef GRF_FORESTTRAINER_H
#define GRF_FORESTTRAINER_H

#include <memory>

#include "prediction/OptimizedPredictionStrategy.h"
#include "relabeling/RelabelingStrategy.h"
#include "splitting/factory/SplittingRuleFactory.h"

#include "tree/Tree.h"
#include "tree/TreeTrainer.h"
#include "forest/Forest.h"
#include "ForestOptions.h"

namespace grf {

class ForestTrainer {
public:
  ForestTrainer(std::unique_ptr<RelabelingStrategy> relabeling_strategy,
                std::unique_ptr<SplittingRuleFactory> splitting_rule_factory,
                std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy);

  Forest train(const Data& data, const ForestOptions& options) const;

private:

  std::vector<std::unique_ptr<Tree>> train_trees(const Data& data,
                                                 const ForestOptions& options) const;

  std::vector<std::unique_ptr<Tree>> train_batch(
      size_t start,
      size_t num_trees,
      const Data& data,
      const ForestOptions& options) const;

  std::unique_ptr<Tree> train_tree(const Data& data,
                                   RandomSampler& sampler,
                                   const ForestOptions& options) const;

  std::vector<std::unique_ptr<Tree>> train_ci_group(const Data& data,
                                                    RandomSampler& sampler,
                                                    const ForestOptions& options) const;

  TreeTrainer tree_trainer;
};

} // namespace grf

#endif //GRF_FORESTTRAINER_H
