
#include <iterator>
#include "BootstrapSampler.h"

#include "Tree.h"
#include "utility.h"

Tree::Tree(const std::vector<std::vector<size_t>> &child_nodeIDs,
           const std::vector<std::vector<size_t>> sampleIDs,
           const std::vector<size_t> &split_varIDs,
           const std::vector<double> &split_values,
           const std::vector<size_t> oob_sampleIDs,
           const std::vector<size_t> inbag_counts) :
    child_nodeIDs(child_nodeIDs),
    sampleIDs(sampleIDs),
    split_varIDs(split_varIDs),
    split_values(split_values),
    oob_sampleIDs(oob_sampleIDs),
    inbag_counts(inbag_counts) {}

Tree::~Tree() {}

std::vector<size_t> Tree::get_terminal_nodeIDs(const Data *prediction_data,
                                               bool oob_prediction) {
  std::vector<size_t> prediction_terminal_nodeIDs;
  prediction_terminal_nodeIDs.resize(prediction_data->getNumRows());

  size_t num_samples_predict;
  if (oob_prediction) {
    num_samples_predict = oob_sampleIDs.size();
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sampleID;
    if (oob_prediction) {
      sampleID = oob_sampleIDs[i];
    } else {
      sampleID = i;
    }
    size_t nodeID = 0;
    while (1) {

      // Break if terminal node
      if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0) {
        break;
      }

      // Move to child
      size_t split_varID = split_varIDs[nodeID];
      double value = prediction_data->get(sampleID, split_varID);
      if (value <= split_values[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[0][nodeID];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[1][nodeID];
      }
    }

    prediction_terminal_nodeIDs[sampleID] = nodeID;
  }
  return prediction_terminal_nodeIDs;
}


