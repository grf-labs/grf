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

#include <unordered_map>
#include <random>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include "TreeClassification.h"
#include "utility.h"
#include "Data.h"

TreeClassification::TreeClassification(std::vector<double>* class_values, std::vector<uint>* response_classIDs) :
    class_values(class_values), response_classIDs(response_classIDs) {
}

TreeClassification::TreeClassification(std::vector<std::vector<size_t>>& child_nodeIDs,
    std::vector<size_t>& split_varIDs, std::vector<double>& split_values, std::vector<double>* class_values,
    std::vector<uint>* response_classIDs) :
    Tree(child_nodeIDs, split_varIDs, split_values), class_values(class_values), response_classIDs(response_classIDs) {
}

TreeClassification::~TreeClassification() {
}

void TreeClassification::initInternal() {
  // TODO
}

void TreeClassification::addPrediction(size_t nodeID, size_t sampleID) {
  predictions[0][sampleID] = split_values[nodeID];
}

double TreeClassification::estimate(size_t nodeID) {

  // Count classes over samples in node and return class with maximum count
  std::unordered_map<double, size_t> class_count;
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    double value = data->get(sampleIDs[nodeID][i], dependent_varID);
    ++class_count[value];
  }
  return (mostFrequentValue(class_count, random_number_generator));
}

void TreeClassification::appendToFileInternal(std::ofstream& file) {
  // Empty on purpose
}

bool TreeClassification::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  // TODO: <= okay? should it not be < ?
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

  // Find best split, stop if no decrease of impurity
  bool stop = findBestSplit(nodeID, possible_split_varIDs);
  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}

void TreeClassification::createEmptyNodeInternal() {
  // Empty on purpose
}

double TreeClassification::computePredictionAccuracyInternal() {

  size_t num_predictions = predictions[0].size();
  size_t num_missclassifications = 0;
  for (size_t i = 0; i < num_predictions; ++i) {
    double predicted_value = predictions[0][i];
    double real_value = data->get(oob_sampleIDs[i], dependent_varID);
    if (predicted_value != real_value) {
      ++num_missclassifications;
    }
  }
  return (1.0 - (double) num_missclassifications / (double) num_predictions);
}

bool TreeClassification::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  size_t num_classes = class_values->size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Class counts
  size_t* class_counts_left = new size_t[num_classes];
  size_t* class_counts = new size_t[num_classes];

  // Reuse class_counts_left
  size_t* class_counts_0 = class_counts_left;
  size_t* class_counts_1 = new size_t[num_classes];

  // Compute overall class counts
  for (size_t i = 0; i < num_classes; ++i) {
    class_counts[i] = 0;
  }
  for (size_t i = 0; i < num_samples_node; ++i) {
    size_t sampleID = sampleIDs[nodeID][i];
    uint sample_classID = (*response_classIDs)[sampleID];
    ++class_counts[sample_classID];
  }

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {

    // Create possible split values
    std::vector<double> possible_split_values;
    data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

    //Try next variable if all equal for this
    if (possible_split_values.size() < 2) {
      continue;
    }

    // If the column consists only of gwa data
    if ((possible_split_values.size() == 2 && possible_split_values[0] == 0 && possible_split_values[1] == 1)
        || (possible_split_values.size() == 3 && possible_split_values[0] == 0 && possible_split_values[1] == 1
            && possible_split_values[2] == 2)) {

      findBestSplitValueGWA(nodeID, varID, num_classes, num_samples_node, class_counts, class_counts_0, class_counts_1,
          best_value, best_varID, best_decrease);

    } else {

      findBestSplitValue(nodeID, varID, possible_split_values, class_counts, class_counts_left, num_classes,
          num_samples_node, best_value, best_varID, best_decrease);

    }
  }

  delete[] class_counts_left;
  delete[] class_counts;
  delete[] class_counts_1;

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Compute gini index for this node and to variable importance if needed
  if (importance_mode == IMP_GINI) {
    addGiniImportance(nodeID, best_varID, best_decrease);
  }
  return false;
}

void TreeClassification::findBestSplitValue(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
    size_t* class_counts, size_t* class_counts_left, size_t num_classes, size_t num_samples_node, double& best_value,
    size_t& best_varID, double& best_decrease) {

  // For all possible split values
  for (auto& split_value : possible_split_values) {

    // Initialize
    size_t n_left = 0;
    for (size_t i = 0; i < num_classes; ++i) {
      class_counts_left[i] = 0;
    }

    // Virtually split at this value. Count overall and for classes.
    for (auto& sampleID : sampleIDs[nodeID]) {
      double value = data->get(sampleID, varID);
      uint sample_classID = (*response_classIDs)[sampleID];
      if (value <= split_value) {
        ++n_left;
        ++class_counts_left[sample_classID];
      }
    }

    // Stop if one child empty
    size_t n_right = num_samples_node - n_left;
    if (n_left == 0 || n_right == 0) {
      continue;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t i = 0; i < num_classes; ++i) {
      size_t class_count_right = class_counts[i] - class_counts_left[i];
      sum_left += class_counts_left[i] * class_counts_left[i];
      sum_right += class_count_right * class_count_right;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = split_value;
      best_varID = varID;
      best_decrease = decrease;
    }
  }

}

void TreeClassification::findBestSplitValueGWA(size_t nodeID, size_t varID, size_t num_classes, size_t num_samples_node,
    size_t* class_counts, size_t* class_counts_0, size_t* class_counts_1, double& best_value, size_t& best_varID,
    double& best_decrease) {

  // Initialize
  for (size_t i = 0; i < num_classes; ++i) {
    class_counts_0[i] = 0;
    class_counts_1[i] = 0;
  }
  size_t n_0 = 0;
  size_t n_1 = 0;

  // Count 0,1,2 per classes
  for (size_t i = 0; i < num_samples_node; ++i) {
    size_t sampleID = sampleIDs[nodeID][i];
    double value = data->get(sampleID, varID);
    uint sample_classID = (*response_classIDs)[sampleID];

    // Add to count
    if (value == 0) {
      ++class_counts_0[sample_classID];
      ++n_0;
    } else if (value == 1) {
      ++class_counts_1[sample_classID];
      ++n_1;
    }
  }

  // Split at 0 or 1
  // Sum of squares
  double sum_left_split_0 = 0;
  double sum_right_split_0 = 0;
  double sum_left_split_1 = 0;
  double sum_right_split_1 = 0;
  for (size_t i = 0; i < num_classes; ++i) {
    size_t class_count_right = class_counts[i] - class_counts_0[i];
    sum_left_split_0 += class_counts_0[i] * class_counts_0[i];
    sum_right_split_0 += class_count_right * class_count_right;

    size_t class_count_left = class_counts_0[i] + class_counts_1[i];
    class_count_right = class_counts[i] - class_count_left;
    sum_left_split_1 += class_count_left * class_count_left;
    sum_right_split_1 += class_count_right * class_count_right;
  }

  // Number of samples left and right per split
  size_t n_left_split_0 = n_0;
  size_t n_right_split_0 = num_samples_node - n_0;
  size_t n_left_split_1 = n_0 + n_1;
  size_t n_right_split_1 = num_samples_node - n_left_split_1;

  // Decrease of impurity
  double decrease_split_0 = sum_left_split_0 / (double) n_left_split_0 + sum_right_split_0 / (double) n_right_split_0;
  double decrease_split_1 = sum_left_split_1 / (double) n_left_split_1 + sum_right_split_1 / (double) n_right_split_1;

  // If better than before, use this
  if (decrease_split_0 > best_decrease) {
    best_value = 0;
    best_varID = varID;
    best_decrease = decrease_split_0;
  }
  if (decrease_split_1 > best_decrease) {
    best_value = 1;
    best_varID = varID;
    best_decrease = decrease_split_1;
  }
}

void TreeClassification::addGiniImportance(size_t nodeID, size_t varID, double decrease) {

  std::vector<size_t> class_counts;
  class_counts.resize(class_values->size(), 0);

  for (auto& sampleID : sampleIDs[nodeID]) {
    uint sample_classID = (*response_classIDs)[sampleID];
    class_counts[sample_classID]++;
  }
  double sum_node = 0;
  for (auto& class_count : class_counts) {
    sum_node += class_count * class_count;
  }
  double best_gini = decrease - sum_node / (double) sampleIDs[nodeID].size();

  // No variable importance for no split variables
  size_t tempvarID = varID;
  for (auto& skip : *no_split_variables) {
    if (varID >= skip) {
      --tempvarID;
    }
  }
  variable_importance[tempvarID] += best_gini;
}

