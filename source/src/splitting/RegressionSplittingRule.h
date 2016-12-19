#ifndef RANGER_REGRESSIONSPLITTINGRULE_H
#define RANGER_REGRESSIONSPLITTINGRULE_H

#include "Tree.h"
#include "SplittingRule.h"
#include <unordered_map>

class RegressionSplittingRule: public SplittingRule {
public:
  RegressionSplittingRule(Data *data);

  ~RegressionSplittingRule();

  bool findBestSplit(size_t nodeID,
                     std::vector<size_t>& possible_split_varIDs,
                     std::unordered_map<size_t, double> responses_by_sampleID,
                     std::vector<std::vector<size_t>> &sampleIDs,
                     std::vector<size_t> &split_varIDs,
                     std::vector<double> &split_values);

private:
  virtual void findBestSplitValueSmallQ(size_t nodeID,
                                        size_t varID,
                                        double sum_node,
                                        size_t num_samples_node,
                                        double& best_value,
                                        size_t& best_varID,
                                        double& best_decrease,
                                        std::unordered_map<size_t, double> responses_by_sampleID,
                                        std::vector<std::vector<size_t>> &sampleIDs);
  virtual void findBestSplitValueLargeQ(size_t nodeID,
                                        size_t varID,
                                        double sum_node,
                                        size_t num_samples_node,
                                        double& best_value,
                                        size_t& best_varID,
                                        double& best_decrease,
                                        std::unordered_map<size_t, double> responses_by_sampleID,
                                        std::vector<std::vector<size_t>> &sampleIDs);

  Data* data;
  size_t* counter;
  double* sums;

  DISALLOW_COPY_AND_ASSIGN(RegressionSplittingRule);
};


#endif //RANGER_REGRESSIONSPLITTINGRULE_H
