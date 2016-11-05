#ifndef RANGER_TREEQUANTILE_H
#define RANGER_TREEQUANTILE_H

#include "ProbabilitySplittingRule.h"
#include "RelabelingStrategy.h"
#include "TreeFactory.h"

class QuantileTreeFactory : public TreeFactory {
public:
  QuantileTreeFactory(RelabelingStrategy *relabeling_strategy, SplittingRule* splitting_rule);

  QuantileTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                      std::vector<double> &split_values,
                      RelabelingStrategy *relabeling_strategy,
                      SplittingRule* splitting_rule,
                      std::vector<std::vector<size_t>> sampleIDs);

  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  std::vector<size_t> get_neighboring_samples(size_t sampleID);

private:
  RelabelingStrategy* relabeling_strategy;
  SplittingRule* splitting_rule;

  DISALLOW_COPY_AND_ASSIGN(QuantileTreeFactory);
};
#endif //RANGER_TREEQUANTILE_H
