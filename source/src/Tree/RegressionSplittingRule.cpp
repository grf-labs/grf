#include "RegressionSplittingRule.h"

RegressionSplittingRule::RegressionSplittingRule(std::unordered_map<size_t, double> responses_by_sampleID,
                                                 Data *data, std::vector<std::vector<size_t>> &sampleIDs,
                                                 std::vector<size_t> &split_varIDs,
                                                 std::vector<double> &split_values) :
    split_varIDs(split_varIDs), split_values(split_values), sampleIDs(sampleIDs) {
  this->responses_by_sampleID = responses_by_sampleID;
  this->data = data;

  size_t max_num_unique_values = data->getMaxNumUniqueValues();
  this->counter = new size_t[max_num_unique_values];
}

RegressionSplittingRule::~RegressionSplittingRule() {
  if (counter != 0) {
    delete[] counter;
  }
  if (sums != 0) {
    delete[] sums;
  }
}

bool RegressionSplittingRule::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  size_t num_samples_node = sampleIDs[nodeID].size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  // Compute sum of responses in node
  double sum_node = 0;
  for (auto& sampleID : sampleIDs[nodeID]) {
    sum_node += responses_by_sampleID[sampleID];
  }

  // For all possible split variables
  for (auto& varID : possible_split_varIDs) {
      // Use faster method for both cases
    double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
    if (q < Q_THRESHOLD) {
      findBestSplitValueSmallQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
    } else {
      findBestSplitValueLargeQ(nodeID, varID, sum_node, num_samples_node, best_value, best_varID, best_decrease);
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

void RegressionSplittingRule::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                                       double &best_value, size_t &best_varID, double &best_decrease) {

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
    sums_right = sums;
    n_right = counter;
    std::fill(sums_right, sums_right + num_splits, 0);
    std::fill(n_right, n_right + num_splits, 0);

  // Sum in right child and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    double response = responses_by_sampleID[sampleID];

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
                                                       double &best_value, size_t &best_varID, double &best_decrease) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill(counter, counter + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);

  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);

    sums[index] += responses_by_sampleID[sampleID];
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