#ifndef RANGER_TreeInstrumental_H
#define RANGER_TreeInstrumental_H

#include "RegressionTreeFactory.h"
#include "RelabelingStrategy.h"

class InstrumentalTreeFactory: public TreeFactory {
public:
  InstrumentalTreeFactory(RelabelingStrategy* relabeling_strategy, SplittingRule* splittingRule);

  InstrumentalTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                          std::vector<double> &split_values,
                          std::vector<std::vector<size_t>> sampleIDs,
                          RelabelingStrategy* relabeling_strategy,
                          SplittingRule* splitting_rule);

  std::vector<size_t> get_neighboring_samples(size_t sampleID);

  bool splitNodeInternal(size_t nodeID, std::vector<size_t> &possible_split_varIDs);

private:
  RelabelingStrategy* relabeling_strategy;
  SplittingRule* splitting_rule;

  DISALLOW_COPY_AND_ASSIGN(InstrumentalTreeFactory);
};
#endif //RANGER_TreeInstrumental_H
