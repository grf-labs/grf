/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 wright@imbs.uni-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <iterator>
#include <Forest/BootstrapSampler.h>

#include "TreeFactory.h"
#include "utility.h"

TreeFactory::TreeFactory(RelabelingStrategy* relabeling_strategy,
                         SplittingRule* splitting_rule) :
        relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule),
        dependent_varID(0), mtry(0), no_split_variables(0), min_node_size(0),
        deterministic_varIDs(0), split_select_varIDs(0), split_select_weights(0), data(0) {}

TreeFactory::TreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs,
                         std::vector<size_t> &split_varIDs,
                         std::vector<double> &split_values,
                         std::vector<std::vector<size_t>> sampleIDs,
                         RelabelingStrategy *relabeling_strategy,
                         SplittingRule *splitting_rule) :
    dependent_varID(0), mtry(0), no_split_variables(
    0), min_node_size(0), deterministic_varIDs(0), split_select_varIDs(0), split_select_weights(0),
    split_varIDs(split_varIDs), split_values(split_values), child_nodeIDs(child_nodeIDs),
    relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule), data(0) {}

TreeFactory::~TreeFactory() {}

void TreeFactory::init(Data* data, uint mtry, size_t dependent_varID, size_t num_samples, uint seed,
    std::vector<size_t>* deterministic_varIDs, std::vector<size_t>* split_select_varIDs,
    std::vector<double>* split_select_weights, uint min_node_size,
    std::vector<size_t>* no_split_variables, bool sample_with_replacement,
    std::vector<double>* case_weights, bool keep_inbag,
    double sample_fraction) {

  this->data = data;
  this->mtry = mtry;
  this->dependent_varID = dependent_varID;

  // Create root node, assign bootstrap sample and oob samples
  child_nodeIDs.push_back(std::vector<size_t>());
  child_nodeIDs.push_back(std::vector<size_t>());
  createEmptyNode();

  this->deterministic_varIDs = deterministic_varIDs;
  this->split_select_varIDs = split_select_varIDs;
  this->split_select_weights = split_select_weights;
  this->min_node_size = min_node_size;
  this->no_split_variables = no_split_variables;

  this->bootstrap_sampler = new BootstrapSampler(num_samples,
      sampleIDs,
      seed,
      sample_with_replacement,
      sample_fraction,
      keep_inbag,
      case_weights);
}

void TreeFactory::grow() {
  bootstrap_sampler->sample();

  // While not all nodes terminal, split next node
  size_t num_open_nodes = 1;
  size_t i = 0;
  while (num_open_nodes > 0) {
    bool is_terminal_node = splitNode(i);
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
}

void TreeFactory::predict(const Data* prediction_data, bool oob_prediction) {

  size_t num_samples_predict;
  if (oob_prediction) {
    num_samples_predict = bootstrap_sampler->getNumSamplesOob();
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

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
}

void TreeFactory::appendToFile(std::ofstream& file) {

// Save general fields
  saveVector2D(child_nodeIDs, file);
  saveVector1D(split_varIDs, file);
  saveVector1D(split_values, file);

  saveVector2D(sampleIDs, file);
}

void TreeFactory::createPossibleSplitVarSubset(std::vector<size_t>& result) {

// Always use deterministic variables
  std::copy(deterministic_varIDs->begin(), deterministic_varIDs->end(), std::inserter(result, result.end()));

// Randomly add non-deterministic variables (according to weights if needed)
  if (split_select_weights->empty()) {
    bootstrap_sampler->drawWithoutReplacementSkip(result, data->getNumCols(), *no_split_variables, mtry);
  } else {
    size_t num_draws = mtry - result.size();
    bootstrap_sampler->drawWithoutReplacementWeighted(result, *split_select_varIDs, num_draws,
        *split_select_weights);
  }
}

bool TreeFactory::splitNode(size_t nodeID) {

// Select random subset of variables to possibly split at
  std::vector<size_t> possible_split_varIDs;
  createPossibleSplitVarSubset(possible_split_varIDs);

// Call subclass method, sets split_varIDs and split_values
  bool stop = splitNodeInternal(nodeID, possible_split_varIDs);
  if (stop) {
    // Terminal node
    return true;
  }

  size_t split_varID = split_varIDs[nodeID];
  double split_value = split_values[nodeID];

// Create child nodes
  size_t left_child_nodeID = sampleIDs.size();
  child_nodeIDs[0][nodeID] = left_child_nodeID;
  createEmptyNode();

  size_t right_child_nodeID = sampleIDs.size();
  child_nodeIDs[1][nodeID] = right_child_nodeID;
  createEmptyNode();

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

void TreeFactory::createEmptyNode() {
  split_varIDs.push_back(0);
  split_values.push_back(0);
  child_nodeIDs[0].push_back(0);
  child_nodeIDs[1].push_back(0);
  sampleIDs.push_back(std::vector<size_t>());
}



bool TreeFactory::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
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

std::vector<size_t> TreeFactory::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}

