#include "TreePredictor.h"

std::vector<size_t> TreePredictor::get_terminal_nodeIDs(std::shared_ptr<Tree> tree,
                                                        Data *prediction_data,
                                                        const std::vector<size_t>& sampleIDs) {
  bool use_subsample = !sampleIDs.empty();
  const std::vector<std::vector<size_t>>& child_nodeIDs = tree->get_child_nodeIDs();
  const std::vector<double>& split_values = tree->get_split_values();
  const std::vector<size_t>& split_varIDs = tree->get_split_varIDs();

  std::vector<size_t> prediction_terminal_nodeIDs;
  prediction_terminal_nodeIDs.resize(prediction_data->getNumRows());

  size_t num_samples_predict;
  if (use_subsample) {
    num_samples_predict = sampleIDs.size();
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sampleID;
    if (use_subsample) {
      sampleID = sampleIDs[i];
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
