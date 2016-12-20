#ifndef GRADIENTFOREST_PROBABILITYSPLITTINGRULE_H
#define GRADIENTFOREST_PROBABILITYSPLITTINGRULE_H

#include <globals.h>
#include <vector>
#include <utility/Data.h>
#include "SplittingRule.h"

class ProbabilitySplittingRule: public SplittingRule {
public:
  ProbabilitySplittingRule(Data* data, size_t num_classes);

  bool findBestSplit(size_t nodeID,
                     std::vector<size_t>& possible_split_varIDs,
                     std::unordered_map<size_t, double> responses_by_sampleID,
                     std::vector<std::vector<size_t>> &sampleIDs,
                     std::vector<size_t> &split_varIDs,
                     std::vector<double> &split_values);

private:
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes, size_t *class_counts,
                                size_t num_samples_node,
                                double &best_value, size_t &best_varID, double &best_decrease,
                                std::unordered_map<size_t, double> response_classIDs,
                                std::vector<std::vector<size_t>> &sampleIDs);

  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes, size_t *class_counts,
                                size_t num_samples_node,
                                double &best_value, size_t &best_varID, double &best_decrease,
                                std::unordered_map<size_t, double> responses_by_sampleID,
                                std::vector<std::vector<size_t>> &sampleIDs);

  Data *data;
  size_t num_classes;

  size_t *counter;
  size_t* counter_per_class;

  DISALLOW_COPY_AND_ASSIGN(ProbabilitySplittingRule);
};

#endif //GRADIENTFOREST_PROBABILITYSPLITTINGRULE_H
