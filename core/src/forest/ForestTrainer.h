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

#include "prediction/OptimizedPredictionStrategy.h"
#include "relabeling/RelabelingStrategy.h"
#include "splitting/factory/SplittingRuleFactory.h"

#include "tree/Tree.h"
#include "tree/TreeTrainer.h"
#include "forest/Forest.h"
#include "ForestOptions.h"

#include <future>
#include <memory>
#include <set>
#include <thread>

class ForestTrainer {
public:
  ForestTrainer(const std::unordered_map<size_t, size_t>& observables,
                std::shared_ptr<RelabelingStrategy> relabeling_strategy,
                std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
                std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy);
  const Forest train(Data* data, const ForestOptions& options) const;

private:

  std::vector<std::shared_ptr<Tree>> train_batch(
      size_t start, size_t num_trees,
      Data* data,
      const Observations& observations,
      const ForestOptions& options) const;

  std::shared_ptr<Tree> train_tree(Data* data,
                                   const Observations& observations,
                                   RandomSampler& sampler,
                                   const ForestOptions& options) const;

  std::vector<std::shared_ptr<Tree>> train_ci_group(Data* data,
                                                    const Observations& observations,
                                                    RandomSampler& sampler,
                                                    const ForestOptions& options) const;

  std::unordered_map<size_t, size_t> observables;
  TreeTrainer tree_trainer;
};


#endif //GRF_FORESTTRAINER_H
