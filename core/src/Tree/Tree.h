#ifndef GRADIENTFOREST_TREE_H_
#define GRADIENTFOREST_TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "globals.h"
#include "Data.h"
#include "sampling/BootstrapSampler.h"
#include "splitting/SplittingRule.h"
#include "relabeling/RelabelingStrategy.h"

class Tree {
public:
  Tree(std::vector<std::vector<size_t>>& child_nodeIDs,
       std::vector<std::vector<size_t>> sampleIDs,
       std::vector<size_t>& split_varIDs,
       std::vector<double>& split_values,
       BootstrapSampler* bootstrap_sampler);

  ~Tree();

  std::vector<size_t> predict(const Data* prediction_data, bool oob_prediction);

  const std::vector<std::vector<size_t> >& getChildNodeIDs() const {
    return child_nodeIDs;
  }

  const std::vector<std::vector<size_t>>& get_sampleIDs() const {
    return sampleIDs;
  }

  const std::vector<double>& getSplitValues() const {
    return split_values;
  }

  const std::vector<size_t>& getSplitVarIDs() const {
    return split_varIDs;
  }

  const std::vector<size_t>& getOobSampleIDs() const {
    return bootstrap_sampler->getOobSampleIDs();
  }

  const std::vector<size_t>& getInbagCounts() const {
    return bootstrap_sampler->getInbagCounts();
  }

private:
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> sampleIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  BootstrapSampler* bootstrap_sampler;

private:
  DISALLOW_COPY_AND_ASSIGN(Tree);
};

#endif /* GRADIENTFOREST_TREE_H_ */
