#ifndef RANGER_TREECAUSAL_H
#define RANGER_TREECAUSAL_H

#include "TreeRegression.h"

class TreeCausal: public TreeRegression {
public:
  TreeCausal(size_t treatment_varID);
  TreeCausal(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
             std::vector<double> &split_values, std::vector<std::vector<size_t>> sampleIDs,
             size_t treatment_varID);
  std::vector<size_t> get_neighboring_samples(size_t sampleID);


  bool splitNodeInternal(size_t nodeID, std::vector<size_t> &possible_split_varIDs);

  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                double &best_value, size_t &best_varID, double &best_decrease,
                                std::unordered_map<size_t, double> &responses_by_sampleID);

  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                double &best_value, size_t &best_varID, double &best_decrease,
                                std::unordered_map<size_t, double> &responses_by_sampleID);

private:
  std::unordered_map<size_t, double> relabelResponses(std::vector<size_t>& nodeSampleIDs);
  double getTerminalNodePrediction(std::vector<size_t>& node_sampleIDs);
  std::unordered_map<bool, double> calculateAverageResponses(std::vector<size_t>& node_sampleIDs);
  bool equalDoubles(double first, double second);

  size_t treatment_varID;
  std::uniform_int_distribution<uint> udist;

  DISALLOW_COPY_AND_ASSIGN(TreeCausal);
};
#endif //RANGER_TREECAUSAL_H
