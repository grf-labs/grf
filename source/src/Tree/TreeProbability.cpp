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
    class_values(class_values), response_classIDs(response_classIDs), counter(0), sums(0) {
}

TreeProbability::TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
    std::vector<std::vector<double>>& terminal_class_counts, std::vector<bool>* is_ordered_variable) :
    Tree(child_nodeIDs, split_varIDs, split_values, is_ordered_variable), class_values(class_values), response_classIDs(
        response_classIDs), terminal_class_counts(terminal_class_counts), counter(0), sums(0) {
}

TreeProbability::~TreeProbability() {
  // Empty on purpose
}

void TreeProbability::initInternal() {
  // Init counters if not in memory efficient mode
  if (!memory_saving_splitting) {
    size_t max_num_unique_values = data->getMaxNumUniqueValues();
    counter = new size_t[max_num_unique_values];
    sums = new double[max_num_unique_values];
  }
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

    // Find best split value, if ordered consider all values as split values, else all 2-partitions
    if ((*is_ordered_variable)[varID]) {

      // Use memory saving method if option set
      if (memory_saving_splitting) {
        findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
      } else {
        // Use faster method for both cases
        double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
        if (q < Q_THRESHOLD) {
          findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
        } else {
          findBestSplitValueLargeQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
        }
      }
    } else {
      findBestSplitValueUnordered(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
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

void TreeProbability::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Remove largest value because no split possible
  possible_split_values.pop_back();

  // Initialize with 0m if not in memory efficient mode, use pre-allocated space
  size_t num_splits = possible_split_values.size();
  double* sums_right;
  size_t* n_right;
  if (memory_saving_splitting) {
    sums_right = new double[num_splits]();
    n_right = new size_t[num_splits]();
  } else {
    sums_right = sums;
    n_right = counter;
    std::fill(sums_right, sums_right + num_splits, 0);
    std::fill(n_right, n_right + num_splits, 0);
  }

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

  if (memory_saving_splitting) {
    delete[] sums_right;
    delete[] n_right;
  }
}

void TreeProbability::findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill(counter, counter + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);

  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);

    sums[index] += data->get(sampleID, dependent_varID);
    ++counter[index];
  }

  size_t n_left = 0;
  double sum_left = 0;

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {

    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];
    sum_left += sums[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    double sum_right = sum_node - sum_left;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = data->getUniqueDataValue(varID, i);
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void TreeProbability::findBestSplitValueUnordered(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
    double& best_value, size_t& best_varID, double& best_decrease) {

  // Create possible split values
  std::vector<double> factor_levels;
  data->getAllValues(factor_levels, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (factor_levels.size() < 2) {
    return;
  }

  // Number of possible splits is 2^num_levels
  size_t num_splits = (1 << factor_levels.size());

  // Compute decrease of impurity for each possible split
  // Split where all left (0) or all right (1) are excluded
  // The second half of numbers is just left/right switched the first half -> Exclude second half
  for (size_t local_splitID = 1; local_splitID < num_splits / 2; ++local_splitID) {

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

