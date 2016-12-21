
#include <iterator>
#include "BootstrapSampler.h"

#include "Tree.h"
#include "utility.h"

Tree::Tree(std::vector<std::vector<size_t>> &child_nodeIDs,
           std::vector<std::vector<size_t>> sampleIDs,
           std::vector<size_t> &split_varIDs,
           std::vector<double> &split_values,
           BootstrapSampler* bootstrap_sampler) :
    child_nodeIDs(child_nodeIDs),
    sampleIDs(sampleIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    bootstrap_sampler(bootstrap_sampler) {}

Tree::~Tree() {}

std::vector<size_t> Tree::predict(const Data* prediction_data, bool oob_prediction) {

  size_t num_samples_predict;
  if (oob_prediction) {
    num_samples_predict = bootstrap_sampler->getNumSamplesOob();
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

  std::vector<size_t> prediction_terminal_nodeIDs;
  prediction_terminal_nodeIDs.resize(num_samples_predict, 0);

// For each sample start in root, drop down the tree and return final value
  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sample_idx;
    if (oob_prediction) {
      sample_idx = bootstrap_sampler->getOobSampleIDs()[i];
    } else {
      sample_idx = i;
    }
    size_t nodeID = 0;
    while (1) {

      // Break if terminal node
      if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0) {
        break;
      }

      // Move to child
      size_t split_varID = split_varIDs[nodeID];
      double value = prediction_data->get(sample_idx, split_varID);
      if (value <= split_values[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

    prediction_terminal_nodeIDs[i] = nodeID;
  }
  return prediction_terminal_nodeIDs;
}


