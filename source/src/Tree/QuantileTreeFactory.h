#ifndef RANGER_TREEQUANTILE_H
#define RANGER_TREEQUANTILE_H

#include "RegressionTreeFactory.h"
#include "ProbabilitySplittingRule.h"
#include "RelabelingStrategy.h"

class QuantileTreeFactory : public TreeFactory {
public:
  QuantileTreeFactory(RelabelingStrategy *relabeling_strategy);

  QuantileTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                      std::vector<double> &split_values,
                      RelabelingStrategy *relabeling_strategy,
                      std::vector<std::vector<size_t>> sampleIDs);

  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  std::vector<size_t> get_neighboring_samples(size_t sampleID);

private:
  ProbabilitySplittingRule* createSplittingRule(std::vector<size_t> &nodeSampleIDs,
                                                std::unordered_map<size_t, double>& relabeledResponses);

  RelabelingStrategy* relabeling_strategy;

  DISALLOW_COPY_AND_ASSIGN(QuantileTreeFactory);
};
#endif //RANGER_TREEQUANTILE_H
