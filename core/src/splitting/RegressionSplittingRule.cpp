/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.

  Authorship: Marvin Wright (wright@imbs.uni-luebeck.de), refactored by
  Julie Tibshirani (jtibs@cs.stanford.edu)
 #-------------------------------------------------------------------------------*/

#include "RegressionSplittingRule.h"

RegressionSplittingRule::RegressionSplittingRule(Data* data) {
  this->data = data;

  size_t max_num_unique_values = data->getMaxNumUniqueValues();
  this->counter = new size_t[max_num_unique_values];
  this->sums = new double[max_num_unique_values];
}

RegressionSplittingRule::~RegressionSplittingRule() {
  if (counter != 0) {
    delete[] counter;
  }
  if (sums != 0) {
    delete[] sums;
  }
}

bool RegressionSplittingRule::findBestSplit(size_t nodeID,
                                            const std::vector<size_t>& possible_split_varIDs,
                                            const std::unordered_map<size_t, double>& labels_by_sampleID,
                                            const std::vector<std::vector<size_t>> &sampleIDs,
                                            std::vector<size_t> &split_varIDs,
                                            std::vector<double> &split_values) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Compute sum of responses in node
  double sum_node = 0;
  for (auto& sampleID : sampleIDs[nodeID]) {
    sum_node += labels_by_sampleID.at(sampleID);
  }

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {
    // Use faster method for both cases
    double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
    if (q < Q_THRESHOLD) {
      findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
                               labels_by_sampleID, sampleIDs);
    } else {
      findBestSplitValueLargeQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease,
                               labels_by_sampleID, sampleIDs);
    }
  }

// Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

// Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;
  return false;
}

void RegressionSplittingRule::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node,
                                                       size_t num_samples_node,
                                                       double &best_value, size_t &best_varID, double &best_decrease,
                                                       const std::unordered_map<size_t, double>& responses_by_sampleID,
                                                       const std::vector<std::vector<size_t>>& sampleIDs) {
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs.at(nodeID), varID);

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
    sums_right = sums;
    n_right = counter;
    std::fill(sums_right, sums_right + num_splits, 0);
    std::fill(n_right, n_right + num_splits, 0);

  // Sum in right child and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    double response = responses_by_sampleID.at(sampleID);

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
}

void RegressionSplittingRule::findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                                       double &best_value, size_t &best_varID, double &best_decrease,
                                                       const std::unordered_map<size_t, double>& responses_by_sampleID,
                                                       const std::vector<std::vector<size_t>> &sampleIDs) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill(counter, counter + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);

  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);

    sums[index] += responses_by_sampleID.at(sampleID);
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