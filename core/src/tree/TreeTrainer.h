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

#ifndef GRADIENTFOREST_TREETRAINER_H
#define GRADIENTFOREST_TREETRAINER_H

#include <memory>

#include "PredictionStrategy.h"
#include "RelabelingStrategy.h"
#include "SplittingRuleFactory.h"
#include "Data.h"
#include "Tree.h"
#include "TreeOptions.h"
#include "BootstrapSampler.h"

class TreeTrainer {
public:
  TreeTrainer(std::shared_ptr<RelabelingStrategy> relabeling_strategy,
              std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
              std::shared_ptr<PredictionStrategy> prediction_strategy,
              const TreeOptions& options);

  std::shared_ptr<Tree> train(Data* data,
                              const Observations& observations,
                              BootstrapSampler& bootstrap_sampler,
                              const std::vector<size_t>& sampleIDs);

private:
  void create_empty_node(std::vector<std::vector<size_t>>& child_nodeIDs,
                         std::vector<std::vector<size_t>>& sampleIDs,
                         std::vector<size_t>& split_varIDs,
                         std::vector<double>& split_values);

  void repopulate_terminal_nodeIDs(std::shared_ptr<Tree> tree,
                                   Data* data,
                                   const std::vector<size_t>& leaf_sampleIDs);

  void create_split_variable_subset(std::vector<size_t>& result,
                                    BootstrapSampler &bootstrap_sampler,
                                    Data *data,
                                    const std::vector<double>& split_select_weights);

  bool split_node(size_t nodeID,
                  std::shared_ptr<SplittingRule> splitting_rule,
                  BootstrapSampler& bootstrap_sampler,
                  Data* data,
                  const Observations& observations,
                  std::vector<std::vector<size_t>>& child_nodeIDs,
                  std::vector<std::vector<size_t>>& sampleIDs,
                  std::vector<size_t>& split_varIDs,
                  std::vector<double>& split_values,
                  const std::vector<double>& split_select_weights);

  bool split_node_internal(size_t nodeID,
                           std::shared_ptr<SplittingRule> splitting_rule,
                           const Observations& observations,
                           const std::vector<size_t>& possible_split_varIDs,
                           std::vector<std::vector<size_t>>& sampleIDs,
                           std::vector<size_t>& split_varIDs,
                           std::vector<double>& split_values);

  std::shared_ptr<RelabelingStrategy> relabeling_strategy;
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory;
  std::shared_ptr<PredictionStrategy> prediction_strategy;

  TreeOptions options;
};

#endif //GRADIENTFOREST_TREETRAINER_H
