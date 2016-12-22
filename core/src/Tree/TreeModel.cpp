
#include "Data.h"
#include "TreeModel.h"

TreeModel::TreeModel(RelabelingStrategy *relabeling_strategy,
                     SplittingRule *splitting_rule,
                     PredictionStrategy *prediction_strategy,
                     size_t dependent_varID,
                     uint mtry,
                     uint min_node_size,
                     std::vector<size_t>* deterministic_varIDs,
                     std::vector<size_t>* split_select_varIDs,
                     std::vector<size_t>* no_split_variables):
    relabeling_strategy(relabeling_strategy),
    splitting_rule(splitting_rule),
    prediction_strategy(prediction_strategy),
    dependent_varID(dependent_varID),
    mtry(mtry),
    min_node_size(min_node_size),
    deterministic_varIDs(deterministic_varIDs),
    split_select_varIDs(split_select_varIDs),
    no_split_variables(no_split_variables) {}


Tree* TreeModel::train(Data* data,
                       BootstrapSampler* bootstrap_sampler,
                       std::unordered_map<std::string, std::vector<double>>* observations,
                       std::vector<double>* split_select_weights) {
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> sampleIDs;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  child_nodeIDs.push_back(std::vector<size_t>());
  child_nodeIDs.push_back(std::vector<size_t>());
  createEmptyNode(child_nodeIDs, sampleIDs, split_varIDs, split_values);

  bootstrap_sampler->sample(sampleIDs);

  // While not all nodes terminal, split next node
  size_t num_open_nodes = 1;
  size_t i = 0;
  while (num_open_nodes > 0) {
    bool is_terminal_node = splitNode(i, bootstrap_sampler,
                                      data, observations, child_nodeIDs,
                                      sampleIDs,
                                      split_varIDs,
                                      split_values,
                                      split_select_weights);
    if (is_terminal_node) {
      --num_open_nodes;
    } else {
      sampleIDs[i].clear();
      ++num_open_nodes;
    }
    ++i;
  }

// Delete sampleID vector to save memory
  //sampleIDs.clear();

  return new Tree(child_nodeIDs,
                  sampleIDs,
                  split_varIDs,
                  split_values,
                  bootstrap_sampler);
}

void TreeModel::createPossibleSplitVarSubset(std::vector<size_t> &result,
                                             BootstrapSampler* bootstrap_sampler,
                                             Data* data,
                                             std::vector<double>* split_select_weights) {

// Always use deterministic variables
  std::copy(deterministic_varIDs->begin(), deterministic_varIDs->end(), std::inserter(result, result.end()));

// Randomly add non-deterministic variables (according to weights if needed)
  if (split_select_weights->empty()) {
    bootstrap_sampler->drawWithoutReplacementSkip(result,
                                                 data->getNumCols(), *no_split_variables, mtry);
  } else {
    size_t num_draws = mtry - result.size();
    bootstrap_sampler->drawWithoutReplacementWeighted(result,
                                                     *split_select_varIDs,
                                                     num_draws,
                                                     *split_select_weights);
  }
}

bool TreeModel::splitNode(size_t nodeID,
                          BootstrapSampler* bootstrap_sampler,
                          Data* data,
                          std::unordered_map<std::string, std::vector<double>>* observations,
                          std::vector<std::vector<size_t>>& child_nodeIDs,
                          std::vector<std::vector<size_t>>& sampleIDs,
                          std::vector<size_t>& split_varIDs,
                          std::vector<double>& split_values,
                          std::vector<double>* split_select_weights) {

// Select random subset of variables to possibly split at
  std::vector<size_t> possible_split_varIDs;
  createPossibleSplitVarSubset(possible_split_varIDs, bootstrap_sampler, data, split_select_weights);

// Call subclass method, sets split_varIDs and split_values
  bool stop = splitNodeInternal(nodeID,
                                data,
                                observations,
                                possible_split_varIDs,
                                sampleIDs,
                                split_varIDs,
                                split_values);
  if (stop) {
    // Terminal node
    return true;
  }

  size_t split_varID = split_varIDs[nodeID];
  double split_value = split_values[nodeID];

// Create child nodes
  size_t left_child_nodeID = sampleIDs.size();
  child_nodeIDs[0][nodeID] = left_child_nodeID;
  createEmptyNode(child_nodeIDs, sampleIDs, split_varIDs, split_values);

  size_t right_child_nodeID = sampleIDs.size();
  child_nodeIDs[1][nodeID] = right_child_nodeID;
  createEmptyNode(child_nodeIDs, sampleIDs, split_varIDs, split_values);

// For each sample in node, assign to left or right child
  // Ordered: left is <= splitval and right is > splitval
  for (auto &sampleID : sampleIDs[nodeID]) {
    if (data->get(sampleID, split_varID) <= split_value) {
      sampleIDs[left_child_nodeID].push_back(sampleID);
    } else {
      sampleIDs[right_child_nodeID].push_back(sampleID);
    }
  }

// No terminal node
  return false;
}

bool TreeModel::splitNodeInternal(size_t nodeID,
                                  Data* data,
                                  std::unordered_map<std::string, std::vector<double>>* observations,
                                  std::vector<size_t>& possible_split_varIDs,
                                  std::vector<std::vector<size_t>>& sampleIDs,
                                  std::vector<size_t>& split_varIDs,
                                  std::vector<double>& split_values) {
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

  std::unordered_map<size_t, double> responses_by_sampleID = relabeling_strategy->relabel_outcomes(
      observations, sampleIDs[nodeID]);

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

void TreeModel::createEmptyNode(std::vector<std::vector<size_t>>& child_nodeIDs,
                                std::vector<std::vector<size_t>>& sampleIDs,
                                std::vector<size_t>& split_varIDs,
                                std::vector<double>& split_values) {
  split_varIDs.push_back(0);
  split_values.push_back(0);
  child_nodeIDs[0].push_back(0);
  child_nodeIDs[1].push_back(0);
  sampleIDs.push_back(std::vector<size_t>());
}