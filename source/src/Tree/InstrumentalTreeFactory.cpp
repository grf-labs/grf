
#include <set>
#include <unordered_map>
#include "utility.h"
#include "InstrumentalTreeFactory.h"
#include "InstrumentalRelabelingStrategy.h"

InstrumentalTreeFactory::InstrumentalTreeFactory(RelabelingStrategy* relabeling_strategy,
                                                 SplittingRule* splitting_rule) :
    relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule) {}

InstrumentalTreeFactory::InstrumentalTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs,
                                                 std::vector<size_t> &split_varIDs,
                                                 std::vector<double> &split_values,
                                                 std::vector<std::vector<size_t>> sampleIDs,
                                                 RelabelingStrategy *relabeling_strategy,
                                                 SplittingRule* splitting_rule) :
    TreeFactory(child_nodeIDs, split_varIDs, split_values),
    relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule) {
  this->sampleIDs = sampleIDs;
}

bool InstrumentalTreeFactory::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    split_values[nodeID] = -1.0;
    return true;
  }

  // Check if node is pure and set split_value to estimate and stop if pure
  bool pure = true;
  double pure_value = 0;
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    double value = data->get(sampleIDs[nodeID][i], dependent_varID);
    if (i != 0 && value != pure_value) {
      pure = false;
      break;
    }
    pure_value = value;
  }

  if (pure) {
    split_values[nodeID] = -1.0;
    return true;
  }

  std::unordered_map<size_t, double> responses_by_sampleID = relabeling_strategy->relabelResponses(
      data, sampleIDs[nodeID]);

  bool stop = responses_by_sampleID.empty() ||
      splitting_rule->findBestSplit(nodeID,
                                    possible_split_varIDs,
                                    responses_by_sampleID,
                                    sampleIDs,
                                    split_varIDs,
                                    split_values);

  if (stop) {
    split_values[nodeID] = -1.0;
    return true;
  }
  return false;
}

std::vector<size_t> InstrumentalTreeFactory::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}

