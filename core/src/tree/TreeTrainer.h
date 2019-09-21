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

namespace grf {

class TreeTrainer {
public:
  TreeTrainer(std::unique_ptr<RelabelingStrategy> relabeling_strategy,
              std::unique_ptr<SplittingRuleFactory> splitting_rule_factory,
              std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy);

  std::unique_ptr<Tree> train(const Data& data,
                              RandomSampler& sampler,
                              const std::vector<size_t>& clusters,
                              const TreeOptions& options) const;

private:
  void create_empty_node(std::vector<std::vector<size_t>>& child_nodes,
                         std::vector<std::vector<size_t>>& samples,
                         std::vector<size_t>& split_vars,
                         std::vector<double>& split_values) const;

  void repopulate_leaf_nodes(const std::unique_ptr<Tree>& tree,
                             const Data& data,
                             const std::vector<size_t>& leaf_samples,
                             const bool honesty_prune_leaves) const;

  void create_split_variable_subset(std::vector<size_t>& result,
                                    RandomSampler& sampler,
                                    const Data& data,
                                    uint mtry) const;

  bool split_node(size_t node,
                  const Data& data,
                  const std::unique_ptr<SplittingRule>& splitting_rule,
                  RandomSampler& sampler,
                  std::vector<std::vector<size_t>>& child_nodes,
                  std::vector<std::vector<size_t>>& samples,
                  std::vector<size_t>& split_vars,
                  std::vector<double>& split_values,
                  std::vector<double>& responses_by_sample,
                  const TreeOptions& tree_options) const;

  bool split_node_internal(size_t node,
                           const Data& data,
                           const std::unique_ptr<SplittingRule>& splitting_rule,
                           const std::vector<size_t>& possible_split_vars,
                           const std::vector<std::vector<size_t>>& samples,
                           std::vector<size_t>& split_vars,
                           std::vector<double>& split_values,
                           std::vector<double>& responses_by_sample,
                           uint min_node_size) const ;

  std::set<size_t> disallowed_split_variables;

  std::unique_ptr<RelabelingStrategy> relabeling_strategy;
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory;
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy;
};

} // namespace grf

#endif //GRF_TREETRAINER_H
