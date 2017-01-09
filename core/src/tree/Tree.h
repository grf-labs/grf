#ifndef GRADIENTFOREST_TREE_H_
#define GRADIENTFOREST_TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "globals.h"
#include "Data.h"
#include "BootstrapSampler.h"
#include "SplittingRule.h"
#include "RelabelingStrategy.h"

class Tree {
public:
  Tree(const std::vector<std::vector<size_t>>& child_nodeIDs,
       const std::vector<std::vector<size_t>> sampleIDs,
       const std::vector<size_t>& split_varIDs,
       const std::vector<double>& split_values,
       const std::vector<size_t> oob_sampleIDs,
       const std::vector<size_t> inbag_counts);

  ~Tree();

  std::vector<size_t> get_terminal_nodeIDs(const Data* prediction_data,
                                           bool oob_prediction);

  std::vector<std::vector<size_t> >& getChildNodeIDs() {
    return child_nodeIDs;
  }

  std::vector<std::vector<size_t>>& get_sampleIDs() {
    return sampleIDs;
  }

  std::vector<double>& getSplitValues() {
    return split_values;
  }

  std::vector<size_t>& getSplitVarIDs() {
    return split_varIDs;
  }

  std::vector<size_t>& getOobSampleIDs() {
    return oob_sampleIDs;
  }

  std::vector<size_t>& get_inbag_counts() {
    return inbag_counts;
  }

private:
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> sampleIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  std::vector<size_t> oob_sampleIDs;
  std::vector<size_t> inbag_counts;
};

#endif /* GRADIENTFOREST_TREE_H_ */
