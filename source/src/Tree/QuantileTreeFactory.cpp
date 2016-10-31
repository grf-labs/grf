
#include <set>
#include "QuantileTreeFactory.h"
#include "utility.h"

QuantileTreeFactory::QuantileTreeFactory(RelabelingStrategy* relabeling_strategy) : relabeling_strategy(relabeling_strategy) {}

QuantileTreeFactory::QuantileTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs,
                                         std::vector<size_t> &split_varIDs,
                                         std::vector<double> &split_values,
                                         RelabelingStrategy *relabeling_strategy,
                                         std::vector<std::vector<size_t>> sampleIDs) :
    TreeFactory(child_nodeIDs, split_varIDs, split_values),
    relabeling_strategy(relabeling_strategy) {
  this->sampleIDs = sampleIDs;
}

bool QuantileTreeFactory::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
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

  std::unordered_map<size_t, double> relabeled_responses = relabeling_strategy->relabelResponses(data, sampleIDs[nodeID]);

  ProbabilitySplittingRule* splittingRule = createSplittingRule(sampleIDs[nodeID], relabeled_responses);
  bool stop = splittingRule->findBestSplit(nodeID, possible_split_varIDs);

  if (stop) {
    split_values[nodeID] = -1.0;
    return true;
  }

  return false;
}

ProbabilitySplittingRule* QuantileTreeFactory::createSplittingRule(std::vector<size_t> &nodeSampleIDs,
                                                                   std::unordered_map<size_t, double>& relabeledResponses) {
  std::unordered_map<size_t, uint> relabeled_classIDs;
  std::set<uint> unique_classIDs;
  for (auto& entry : relabeledResponses) {
    uint classID = (uint) round(entry.second);

    relabeled_classIDs[entry.first] = classID;
    unique_classIDs.insert(classID);
  }
  size_t num_classes = unique_classIDs.size();

  return new ProbabilitySplittingRule(num_classes, relabeled_classIDs,
      data, sampleIDs, split_varIDs, split_values);
}

std::vector<size_t> QuantileTreeFactory::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}
