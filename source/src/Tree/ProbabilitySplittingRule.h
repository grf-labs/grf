
#ifndef RANGER_PROBABILITYSPLITTINGRULE_H
#define RANGER_PROBABILITYSPLITTINGRULE_H

#include <globals.h>
#include <vector>
#include <utility/Data.h>

class ProbabilitySplittingRule {
public:
  ProbabilitySplittingRule(size_t num_classes,
                           std::vector<uint>* response_classIDs,
                           Data* data,
                           std::vector<std::vector<size_t>> &sampleIDs,
                           std::vector<size_t> &split_varIDs,
                           std::vector<double> &split_values);

  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts, size_t num_samples_node,
                                double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts, size_t num_samples_node,
                                double& best_value, size_t& best_varID, double& best_decrease);

private:
  size_t num_classes;
  std::vector<uint>* response_classIDs;

  size_t* counter;
  size_t* counter_per_class;

  bool memory_saving_splitting;

  Data* data;
  std::vector<size_t>& split_varIDs;

  // Value to split at for each node, for now only binary split
  // For terminal nodes the prediction value is saved here
  std::vector<double>& split_values;

  // For each node a vector with IDs of samples in node
  std::vector<std::vector<size_t>>& sampleIDs;
};

#endif //RANGER_PROBABILITYSPLITTINGRULE_H
