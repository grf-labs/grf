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
       const std::vector<std::vector<size_t>>& terminal_nodeIDs,
       const std::vector<size_t>& split_varIDs,
       const std::vector<double>& split_values,
       const std::vector<size_t>& oob_sampleIDs,
       const std::vector<size_t>& inbag_counts);

  ~Tree();

  std::vector<std::vector<size_t> >& get_child_nodeIDs() {
    return child_nodeIDs;
  }

  std::vector<std::vector<size_t>>& get_terminal_nodeIDs() {
    return terminal_nodeIDs;
  }

  std::vector<double>& get_split_values() {
    return split_values;
  }

  std::vector<size_t>& get_split_varIDs() {
    return split_varIDs;
  }

  std::vector<size_t>& get_oob_sampleIDs() {
    return oob_sampleIDs;
  }

  std::vector<size_t>& get_inbag_counts() {
    return inbag_counts;
  }

  void set_terminal_nodeIDs(const std::vector<std::vector<size_t>>& terminal_nodeIDs) {
    this->terminal_nodeIDs = terminal_nodeIDs;
  }

private:
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> terminal_nodeIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  std::vector<size_t> oob_sampleIDs;
  std::vector<size_t> inbag_counts;
};

#endif /* GRADIENTFOREST_TREE_H_ */
