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

#include "Tree.h"
#include "utility.h"

Tree::Tree() :
    dependent_varID(0), mtry(0), num_samples(0), num_samples_oob(0), no_split_variables(0), min_node_size(0), deterministic_varIDs(
        0), split_select_varIDs(0), split_select_weights(0), oob_sampleIDs(0), data(0), importance_mode(
        DEFAULT_IMPORTANCE_MODE), sample_with_replacement(true), splitrule(DEFAULT_SPLITRULE) {
}

Tree::Tree(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values) :
    dependent_varID(0), mtry(0), num_samples(0), num_samples_oob(0), no_split_variables(0), min_node_size(0), deterministic_varIDs(
        0), split_select_varIDs(0), split_select_weights(0), split_varIDs(split_varIDs), split_values(split_values), child_nodeIDs(
        child_nodeIDs), oob_sampleIDs(0), data(0), importance_mode(DEFAULT_IMPORTANCE_MODE), sample_with_replacement(
        true), splitrule(DEFAULT_SPLITRULE) {
}

Tree::~Tree() {
}

void Tree::init(Data* data, uint mtry, size_t dependent_varID, size_t num_samples, uint seed,
    std::vector<size_t>* deterministic_varIDs, std::vector<size_t>* split_select_varIDs,
    std::vector<double>* split_select_weights, ImportanceMode importance_mode, uint min_node_size,
    std::vector<size_t>* no_split_variables, bool sample_with_replacement, uint splitrule) {

  this->data = data;
  this->mtry = mtry;
  this->dependent_varID = dependent_varID;
  this->num_samples = num_samples;

  // Create root node, assign bootstrap sample and oob samples
  createEmptyNode();

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  this->deterministic_varIDs = deterministic_varIDs;
  this->split_select_varIDs = split_select_varIDs;
  this->split_select_weights = split_select_weights;
  this->importance_mode = importance_mode;
  this->min_node_size = min_node_size;
  this->no_split_variables = no_split_variables;
  this->sample_with_replacement = sample_with_replacement;
  this->splitrule = splitrule;

  // Initialize with variable importance with 0.
  if (importance_mode == IMP_GINI) {
    variable_importance.resize(data->getNumCols() - no_split_variables->size());
  }

  initInternal();
}

void Tree::grow() {

  if (sample_with_replacement) {
    bootstrap();
  } else {
    bootstrapWithoutReplacement();
  }

  // Call recursive split function on root node
  splitNode(0);

  // Delete sampleID vector to save memory
  sampleIDs.clear();
  cleanUpInternal();
}

void Tree::predict(const Data* prediction_data, bool oob_prediction) {

  size_t num_samples_predict;
  if (oob_prediction) {
    num_samples_predict = num_samples_oob;
  } else {
    num_samples_predict = prediction_data->getNumRows();
  }

  reservePredictionMemory(num_samples_predict);

  // For each sample start in root, drop down the tree and return final value
  for (size_t i = 0; i < num_samples_predict; ++i) {
    size_t sample_idx;
    if (oob_prediction) {
      sample_idx = oob_sampleIDs[i];
    } else {
      sample_idx = i;
    }
    size_t nodeID = 0;
    while (1) {

      // Break if terminal node
      if (child_nodeIDs[nodeID].empty()) {
        break;
      }

      // Move to child
      double value = prediction_data->get(sample_idx, split_varIDs[nodeID]);
      if (value <= split_values[nodeID]) {
        // Move to left child
        nodeID = child_nodeIDs[nodeID][0];
      } else {
        // Move to right child
        nodeID = child_nodeIDs[nodeID][1];
      }
    }

    addPrediction(nodeID, i);
  }
}

void Tree::computePermutationImportance() {

  size_t num_independent_variables = data->getNumCols() - no_split_variables->size();
  variable_importance.clear();
  variable_importance.reserve(num_independent_variables);

  // Compute normal prediction accuracy for each tree. Predictions already computed..
  double accuracy_normal = computePredictionAccuracyInternal();

  predictions.clear();
  reservePredictionMemory(num_samples_oob);

  // Reserve space for permutations, initialize with oob_sampleIDs
  std::vector<size_t> permutations(oob_sampleIDs);

  // Randomly permute for all independent variables
  for (size_t i = 0; i < num_independent_variables; ++i) {

    // Skip no split variables
    size_t varID = i;
    for (auto& skip : *no_split_variables) {
      if (varID >= skip) {
        ++varID;
      }
    }

    // Permute and compute prediction accuracy again for this permutation and save difference
    permuteAndPredictOobSamples(varID, permutations);
    double accuracy_permuted = computePredictionAccuracyInternal();
    variable_importance.push_back(accuracy_normal - accuracy_permuted);
  }
}

void Tree::appendToFile(std::ofstream& file) {

  // Save general fields
  saveVector2D(child_nodeIDs, file);
  saveVector1D(split_varIDs, file);
  saveVector1D(split_values, file);

  // Call special functions for subclasses to save special fields.
  appendToFileInternal(file);
}

void Tree::createPossibleSplitVarSubset(std::vector<size_t>& result) {

  // Always use deterministic variables
  std::copy(deterministic_varIDs->begin(), deterministic_varIDs->end(), std::inserter(result, result.end()));

  // Randomly add non-deterministic variables (according to weights if needed)
  if (split_select_weights->empty()) {
    drawWithoutReplacementSkip(result, random_number_generator, data->getNumCols(), *no_split_variables, mtry);
  } else {
    size_t num_draws = mtry - result.size();
    drawWithoutReplacementWeighted(result, random_number_generator, *split_select_varIDs, num_draws,
        *split_select_weights);
  }
}

void Tree::splitNode(size_t nodeID) {

  // Select random subset of variables to possibly split at
  std::vector<size_t> possible_split_varIDs;
  createPossibleSplitVarSubset(possible_split_varIDs);

  // Call subclass method, sets split_varIDs and split_values
  bool stop = splitNodeInternal(nodeID, possible_split_varIDs);
  if (stop) {
    return;
  }

  size_t split_varID = split_varIDs[nodeID];
  double split_value = split_values[nodeID];

  // Create child nodes
  size_t left_child_nodeID = sampleIDs.size();
  child_nodeIDs[nodeID].push_back(left_child_nodeID);
  createEmptyNode();

  size_t right_child_nodeID = sampleIDs.size();
  child_nodeIDs[nodeID].push_back(right_child_nodeID);
  createEmptyNode();

  // For each sample in node, assign to left (<= split val) or right (> split val) child
  for (auto& sampleID : sampleIDs[nodeID]) {
    if (data->get(sampleID, split_varID) <= split_value) {
      sampleIDs[left_child_nodeID].push_back(sampleID);
    } else {
      sampleIDs[right_child_nodeID].push_back(sampleID);
    }
  }

  // Recursively call split node on child nodes
  for (size_t i = 0; i < child_nodeIDs[nodeID].size(); ++i) {
    splitNode(child_nodeIDs[nodeID][i]);
  }
}

void Tree::createEmptyNode() {
  split_varIDs.push_back(0);
  split_values.push_back(0);
  child_nodeIDs.push_back(std::vector<size_t>());
  sampleIDs.push_back(std::vector<size_t>());

  createEmptyNodeInternal();
}

size_t Tree::dropDownSamplePermuted(size_t permuted_varID, size_t sampleID, size_t permuted_sampleID) {

  // Start in root and drop down
  size_t nodeID = 0;
  while (!child_nodeIDs[nodeID].empty()) {

    // Permute if variable is permutation variable
    size_t split_varID = split_varIDs[nodeID];
    size_t sampleID_final = sampleID;
    if (split_varID == permuted_varID) {
      sampleID_final = permuted_sampleID;
    }

    // Move to child
    if (data->get(sampleID_final, split_varID) <= split_values[nodeID]) {
      // Move to left child
      nodeID = child_nodeIDs[nodeID][0];
    } else {
      // Move to right child
      nodeID = child_nodeIDs[nodeID][1];
    }
  }
  return nodeID;
}

void Tree::permuteAndPredictOobSamples(size_t permuted_varID, std::vector<size_t>& permutations) {

  // Permute OOB sample
  //std::vector<size_t> permutations(oob_sampleIDs);
  std::shuffle(permutations.begin(), permutations.end(), random_number_generator);

  // For each sample, drop down the tree and add prediction
  for (size_t i = 0; i < num_samples_oob; ++i) {
    size_t nodeID = dropDownSamplePermuted(permuted_varID, oob_sampleIDs[i], permutations[i]);
    addPrediction(nodeID, i);
  }
}

void Tree::bootstrap() {

  // Reserve space (37% percent are outbag on average, reserve a little more)
  sampleIDs[0].reserve(num_samples);
  oob_sampleIDs.reserve(num_samples * 0.4);

  std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

  // Start with all samples OOB
  std::vector<bool> is_oob;
  is_oob.resize(num_samples, true);

  // Draw num_samples samples with replacement (n out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples; ++s) {
    size_t draw = unif_dist(random_number_generator);
    sampleIDs[0].push_back(draw);
    is_oob[draw] = false;
  }

  // Save OOB samples
  for (size_t s = 0; s < is_oob.size(); ++s) {
    if (is_oob[s]) {
      oob_sampleIDs.push_back(s);
    }
  }
  num_samples_oob = oob_sampleIDs.size();
}

void Tree::bootstrapWithoutReplacement() {

  // As in sampling with replacement use 63.21% of the samples
  size_t num_samples_inbag = (size_t) num_samples * 0.6321;
  shuffleAndSplit(sampleIDs[0], oob_sampleIDs, num_samples, num_samples_inbag, random_number_generator);
  num_samples_oob = oob_sampleIDs.size();
}

