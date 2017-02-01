#ifndef GRADIENTFOREST_TREE_H_
#define GRADIENTFOREST_TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "globals.h"
#include "Data.h"
#include "BootstrapSampler.h"
#include "PredictionValues.h"
#include "SplittingRule.h"

class Tree {
public:
  Tree(const std::vector<std::vector<size_t>>& child_nodeIDs,
       const std::vector<std::vector<size_t>>& leaf_nodeIDs,
       const std::vector<size_t>& split_varIDs,
       const std::vector<double>& split_values,
       const std::vector<size_t>& oob_sampleIDs,
       const PredictionValues& prediction_values);

  ~Tree();

  const std::vector<std::vector<size_t> >& get_child_nodeIDs() {
    return child_nodeIDs;
  }

  const std::vector<std::vector<size_t>>& get_leaf_nodeIDs() {
    return leaf_nodeIDs;
  }

  const std::vector<double>& get_split_values() {
    return split_values;
  }

  const std::vector<size_t>& get_split_varIDs() {
    return split_varIDs;
  }

  const std::vector<size_t>& get_oob_sampleIDs() {
    return oob_sampleIDs;
  }

  const PredictionValues& get_prediction_values() {
    return prediction_values;
  }

  void set_leaf_nodeIDs(const std::vector<std::vector<size_t>>& leaf_nodeIDs) {
    this->leaf_nodeIDs = leaf_nodeIDs;
  }

  void set_oob_sampleIDs(std::vector<size_t>& oob_sampleIDs) {
    this->oob_sampleIDs = oob_sampleIDs;
  }

  void set_prediction_values(const PredictionValues& prediction_values) {
    this->prediction_values = prediction_values;
  }

private:
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> leaf_nodeIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  std::vector<size_t> oob_sampleIDs;

  PredictionValues prediction_values;
};

#endif /* GRADIENTFOREST_TREE_H_ */
