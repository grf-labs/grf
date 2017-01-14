#ifndef GRADIENTFOREST_TREETRAINER_H
#define GRADIENTFOREST_TREETRAINER_H

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
              TreeOptions options);

  std::shared_ptr<Tree> train(Data *data,
                              BootstrapSampler *bootstrap_sampler,
                              const Observations& observations);

private:
  void createEmptyNode(std::vector <std::vector<size_t>> &child_nodeIDs,
                       std::vector <std::vector<size_t>> &sampleIDs,
                       std::vector <size_t> &split_varIDs,
                       std::vector<double> &split_values);

  void repopulate_terminal_nodeIDs(std::shared_ptr<Tree> tree,
                                   Data* data,
                                   std::vector<size_t> leaf_sampleIDs);

  void createPossibleSplitVarSubset(std::vector<size_t> &result,
                                    BootstrapSampler *bootstrap_sampler,
                                    Data *data,
                                    const std::vector<double> &split_select_weights);

  bool splitNode(size_t nodeID,
                 std::shared_ptr<SplittingRule> splitting_rule,
                 BootstrapSampler *bootstrap_sampler,
                 Data *data,
                 const Observations& observations,
                 std::vector <std::vector<size_t>> &child_nodeIDs,
                 std::vector <std::vector<size_t>> &sampleIDs,
                 std::vector <size_t> &split_varIDs,
                 std::vector<double> &split_values,
                 const std::vector<double> &split_select_weights);

  bool splitNodeInternal(size_t nodeID,
                         std::shared_ptr<SplittingRule> splitting_rule,
                         const Observations& observations,
                         const std::vector <size_t> &possible_split_varIDs,
                         std::vector <std::vector<size_t>> &sampleIDs,
                         std::vector <size_t> &split_varIDs,
                         std::vector<double> &split_values);

  std::shared_ptr<RelabelingStrategy> relabeling_strategy;
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory;
  TreeOptions options;
};

#endif //GRADIENTFOREST_TREETRAINER_H
