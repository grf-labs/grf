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

#ifndef GRF_TREETRAINER_H
#define GRF_TREETRAINER_H

#include <memory>

#include "commons/DefaultData.h"
#include "prediction/OptimizedPredictionStrategy.h"
#include "relabeling/RelabelingStrategy.h"
#include "sampling/RandomSampler.h"
#include "splitting/factory/SplittingRuleFactory.h"
#include "tree/Tree.h"
#include "tree/TreeOptions.h"

class TreeTrainer {
public:
  TreeTrainer(std::shared_ptr<RelabelingStrategy> relabeling_strategy,
              std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
              std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy,
              const TreeOptions& options);

  std::shared_ptr<Tree> train(Data* data,
                              const Observations& observations,
                              RandomSampler& sampler,
                              const std::vector<size_t>& samples);

private:
  void create_empty_node(std::vector<std::vector<size_t>>& child_nodes,
                         std::vector<std::vector<size_t>>& samples,
                         std::vector<size_t>& split_vars,
                         std::vector<double>& split_values);

  void repopulate_leaf_nodes(std::shared_ptr<Tree> tree,
                             Data *data,
                             const std::vector<size_t> &leaf_samples);

  void create_split_variable_subset(std::vector<size_t>& result,
                                    RandomSampler& sampler,
                                    Data* data,
                                    const std::vector<double>& split_select_weights);

  bool split_node(size_t node,
                  std::shared_ptr<SplittingRule> splitting_rule,
                  RandomSampler& sampler,
                  Data* data,
                  const Observations& observations,
                  std::vector<std::vector<size_t>>& child_nodes,
                  std::vector<std::vector<size_t>>& samples,
                  std::vector<size_t>& split_vars,
                  std::vector<double>& split_values,
                  const std::vector<double>& split_select_weights);

  bool split_node_internal(size_t node,
                           std::shared_ptr<SplittingRule> splitting_rule,
                           const Observations& observations,
                           const std::vector<size_t>& possible_split_vars,
                           std::vector<std::vector<size_t>>& samples,
                           std::vector<size_t>& split_vars,
                           std::vector<double>& split_values);

  std::shared_ptr<RelabelingStrategy> relabeling_strategy;
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory;
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy;

  TreeOptions options;
};

#endif //GRF_TREETRAINER_H
