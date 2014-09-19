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

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include "TreeRegression.h"
#include "Data.h"

TreeRegression::TreeRegression() {
}

TreeRegression::TreeRegression(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values) :
    Tree(child_nodeIDs, split_varIDs, split_values) {
}

TreeRegression::~TreeRegression() {
}

void TreeRegression::initInternal() {
  // Empty on purpose
}

void TreeRegression::addPrediction(size_t nodeID, size_t sampleID) {
  predictions[0][sampleID] = split_values[nodeID];
}

double TreeRegression::estimate(size_t nodeID) {

  // Mean of responses of samples in node
  double sum_responses_in_node = 0;
  size_t num_samples_in_node = sampleIDs[nodeID].size();
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    sum_responses_in_node += data->get(sampleIDs[nodeID][i], dependent_varID);
  }
  return (sum_responses_in_node / (double) num_samples_in_node);
}

void TreeRegression::appendToFileInternal(std::ofstream& file) {
  // Empty on purpose
}

bool TreeRegression::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() < min_node_size) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  // Find best split, stop if no decrease of impurity
  bool stop = findBestSplit(nodeID, possible_split_varIDs);
  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}

void TreeRegression::createEmptyNodeInternal() {
  // Empty on purpose
}

double TreeRegression::computePredictionAccuracyInternal() {
  double sum_of_squares = 0;
  for (size_t i = 0; i < predictions[0].size(); ++i) {
    double predicted_value = predictions[0][i];
    double real_value = data->get(oob_sampleIDs[i], dependent_varID);
    if (predicted_value != real_value) {
      sum_of_squares += (predicted_value - real_value) * (predicted_value - real_value);
    }
  }
  return (1.0 - sum_of_squares / (double) predictions[0].size());
}

bool TreeRegression::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  //size_t num_samples = sampleIDs[nodeID].size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {

    // Create possible split values
    std::vector<double> possible_split_values;
    data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

    // Try next variable if all equal for this
    if (possible_split_values.size() < 2) {
      continue;
    }

    // For all possible split values
    for (auto& split_value : possible_split_values) {

      // Virtually split at this value. Count and sum overall and for classes.
      size_t n_left = 0;
      size_t n_right = 0;
      double sum_left = 0;
      double sum_right = 0;
      for (auto& sampleID : sampleIDs[nodeID]) {
        double response = data->get(sampleID, dependent_varID);
        if (data->get(sampleID, varID) <= split_value) {
          n_left++;
          sum_left += response;
        } else {
          n_right++;
          sum_right += response;
        }
      }

      // Stop if one child empty
      if (n_left == 0 || n_right == 0) {
        continue;
      }

      // Decrease of impurity = variance reduction = MSE
      double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

      // If better than before, use this
      if (decrease > best_decrease) {
        best_value = split_value;
        best_varID = varID;
        best_decrease = decrease;
      }
    }
  }

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Compute decrease of impurity for this node and add to variable importance if needed
  if (importance_mode == IMP_GINI) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }
  return false;
}

void TreeRegression::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {

  double sum_node = 0;
  for (auto& sampleID : sampleIDs[nodeID]) {
    sum_node += data->get(sampleID, dependent_varID);
  }
  double best_decrease = decrease - sum_node * sum_node / (double) sampleIDs[nodeID].size();

  // No variable importance for no split variables
  size_t tempvarID = varID;
  for (auto& skip : *no_split_variables) {
    if (varID >= skip) {
      --tempvarID;
    }
  }
  variable_importance[tempvarID] += best_decrease;
}

