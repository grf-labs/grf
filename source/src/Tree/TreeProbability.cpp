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

#include "TreeProbability.h"
#include "utility.h"
#include "Data.h"

TreeProbability::TreeProbability(std::vector<double>* class_values, std::vector<uint>* response_classIDs) :
    class_values(class_values), response_classIDs(response_classIDs) {
}

TreeProbability::TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
    std::vector<std::vector<double>>& terminal_class_counts) :
    Tree(child_nodeIDs, split_varIDs, split_values), class_values(class_values), response_classIDs(response_classIDs), terminal_class_counts(
        terminal_class_counts) {
}

TreeProbability::~TreeProbability() {
}

void TreeProbability::initInternal() {
  // Empty on purpose
}

void TreeProbability::addToTerminalNodes(size_t nodeID) {

  size_t num_samples_in_node = sampleIDs[nodeID].size();
  terminal_class_counts[nodeID].resize(class_values->size(), 0);

  // Compute counts
  for (size_t i = 0; i < num_samples_in_node; ++i) {
    size_t node_sampleID = sampleIDs[nodeID][i];
    size_t classID = (*response_classIDs)[node_sampleID];
    ++terminal_class_counts[nodeID][classID];
  }

  // Compute fractions
  for (size_t i = 0; i < terminal_class_counts[nodeID].size(); ++i) {
    terminal_class_counts[nodeID][i] /= num_samples_in_node;
  }
}

void TreeProbability::appendToFileInternal(std::ofstream& file) {

  // Add Terminal node class counts
  // Convert to vector without empty elements and save
  std::vector<size_t> terminal_nodes;
  std::vector<std::vector<double>> terminal_class_counts_vector;
  for (size_t i = 0; i < terminal_class_counts.size(); ++i) {
    if (!terminal_class_counts[i].empty()) {
      terminal_nodes.push_back(i);
      terminal_class_counts_vector.push_back(terminal_class_counts[i]);
    }
  }
  saveVector1D(terminal_nodes, file);
  saveVector2D(terminal_class_counts_vector, file);
}

bool TreeProbability::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    addToTerminalNodes(nodeID);
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
    addToTerminalNodes(nodeID);
    return true;
  }

  // Find best split, stop if no decrease of impurity
  bool stop = findBestSplit(nodeID, possible_split_varIDs);
  if (stop) {
    addToTerminalNodes(nodeID);
    return true;
  }

  return false;
}

void TreeProbability::createEmptyNodeInternal() {
  terminal_class_counts.push_back(std::vector<double>());
}

double TreeProbability::computePredictionAccuracyInternal() {

  size_t num_predictions = prediction_terminal_nodeIDs.size();
  double sum_of_squares = 0;
  for (size_t i = 0; i < num_predictions; ++i) {
    size_t sampleID = oob_sampleIDs[i];
    size_t real_classID = (*response_classIDs)[sampleID];
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    double predicted_value = terminal_class_counts[terminal_nodeID][real_classID];
    sum_of_squares += (1 - predicted_value) * (1 - predicted_value);
  }
  return (1.0 - sum_of_squares / (double) num_predictions);
}

bool TreeProbability::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Compute sum of responses in node
  double sum_node = 0;
  for (auto& sampleID : sampleIDs[nodeID]) {
    sum_node += data->get(sampleID, dependent_varID);
  }

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {

    // Create possible split values
    std::vector<double> all_values;
    data->getAllValues(all_values, sampleIDs[nodeID], varID);

    // Try next variable if all equal for this
    if (all_values.size() < 2) {
      continue;
    }

    // Find best split value, if ordered consider all values as split values, else all 2-partitions
    if ((*is_ordered_variable)[varID]) {

      // Remove largest value because no split possible
      all_values.pop_back();

      findBestSplitValue(nodeID, varID, all_values, sum_node, num_samples_node, best_value, best_varID, best_decrease);
    } else {
      findBestSplitValueUnordered(nodeID, varID, all_values, sum_node, num_samples_node, best_value, best_varID,
          best_decrease);
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

void TreeProbability::findBestSplitValue(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
    double sum_node, size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  size_t num_splits = possible_split_values.size();

  // Initialize with 0
  double* sums_right = new double[num_splits]();
  size_t* n_right = new size_t[num_splits]();

  // Sum in right child and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    double response = data->get(sampleID, dependent_varID);

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        sums_right[i] += response;
      } else {
        break;
      }
    }
  }

  // Compute decrease of impurity for each possible split
  for (size_t i = 0; i < num_splits; ++i) {

    // Stop if one child empty
    size_t n_left = num_samples_node - n_right[i];
    if (n_left == 0 || n_right[i] == 0) {
      continue;
    }

    double sum_right = sums_right[i];
    double sum_left = sum_node - sum_right;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right[i];

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = possible_split_values[i];
      best_varID = varID;
      best_decrease = decrease;
    }
  }

  delete[] sums_right;
  delete[] n_right;

}

void TreeProbability::findBestSplitValueUnordered(size_t nodeID, size_t varID, std::vector<double>& factor_levels,
    double sum_node, size_t num_samples_node, double& best_value, size_t& best_varID, double& best_decrease) {

  // Number of possible splits is 2^num_levels
  size_t num_splits = (1 << factor_levels.size());

  // TODO: Remove
//  std::cout << std::endl << "Number of factors: " << factor_levels.size() << std::endl;

  // Compute decrease of impurity for each possible split
  // Split where all left (0) or all right (1) are excluded
  // The second half of numbers is just left/right switched the first half -> Exclude second half
  for (size_t local_splitID = 1; local_splitID < num_splits / 2; ++local_splitID) {

    // TODO: Remove
//    std::cout << "Local split id: " << local_splitID << std::endl;

    // Compute overall splitID by shifting local factorIDs to global positions
    size_t splitID = 0;
    for (size_t j = 0; j < factor_levels.size(); ++j) {
      if ((local_splitID & (1 << j))) {
        double level = factor_levels[j];
        size_t factorID = floor(level) - 1;
        splitID = splitID | (1 << factorID);
      }
    }

    // Initialize
    double sum_right = 0;
    size_t n_right = 0;

    // Sum in right child
    for (auto& sampleID : sampleIDs[nodeID]) {
      double response = data->get(sampleID, dependent_varID);
      double value = data->get(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count
      // In right child, if bitwise splitID at position factorID is 1
      if ((splitID & (1 << factorID))) {
        ++n_right;
        sum_right += response;
      }
    }
    size_t n_left = num_samples_node - n_right;

    // Sum of squares
    double sum_left = sum_node - sum_right;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

    // TODO: Remove
//    std::cout << "n_left: " << n_left << ", n_right: " << n_right << ", decrease: " << decrease << std::endl;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = splitID;
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void TreeProbability::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {

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

