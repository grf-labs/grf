
#include <set>
#include "TreeQuantile.h"
#include "utility.h"
#include "ProbabilitySplittingRule.h"
#include "QuantileRelabelingStrategy.h"


TreeQuantile::TreeQuantile(std::vector<double>* quantiles) :
    quantiles(quantiles), udist(std::uniform_int_distribution<uint>()) {}

TreeQuantile::TreeQuantile(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                           std::vector<double> &split_values,
                           std::vector<double> *quantiles,
                           std::vector<std::vector<size_t>> sampleIDs) :
    TreeRegression(child_nodeIDs, split_varIDs, split_values),
    quantiles(quantiles),
    udist(std::uniform_int_distribution<uint>()) {
  this->sampleIDs = sampleIDs;
}

bool TreeQuantile::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    split_values[nodeID] = estimate(nodeID);
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
    split_values[nodeID] = pure_value;
    return true;
  }

  std::vector<double> responses;
  for (auto& sampleID : sampleIDs[nodeID]) {
    responses.push_back(data->get(sampleID, dependent_varID));
  }

  QuantileRelabelingStrategy* relabelingStrategy = new QuantileRelabelingStrategy(quantiles);
  std::vector<uint>* relabeled_responses = relabelingStrategy->relabelResponses(&responses);

  ProbabilitySplittingRule* splittingRule = createSplittingRule(sampleIDs[nodeID], relabeled_responses);
  bool stop = splittingRule->findBestSplit(nodeID, possible_split_varIDs);

  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}

ProbabilitySplittingRule* TreeQuantile::createSplittingRule(std::vector<size_t> &nodeSampleIDs,
                                                            std::vector<uint> *relabeledResponses) {
  std::set<uint>* unique_classIDs = new std::set<uint>(relabeledResponses->begin(),
                                                       relabeledResponses->end());
  size_t num_classes = unique_classIDs->size();

  std::vector<uint>* response_classIDs = new std::vector<uint>(num_samples);
  for (size_t i = 0; i < nodeSampleIDs.size(); ++i) {
    size_t sampleID = nodeSampleIDs[i];
    (*response_classIDs)[sampleID] = (*relabeledResponses)[i];
  }

  return new ProbabilitySplittingRule(num_classes, response_classIDs,
      data, sampleIDs, split_varIDs, split_values);
}

std::vector<size_t> TreeQuantile::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}
