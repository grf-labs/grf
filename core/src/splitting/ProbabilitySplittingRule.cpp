#include <unordered_map>
#include "ProbabilitySplittingRule.h"

ProbabilitySplittingRule::ProbabilitySplittingRule(Data* data, size_t num_classes) {
  this->data = data;
  this->num_classes = num_classes;

  size_t max_num_unique_values = data->getMaxNumUniqueValues();
  this->counter = new size_t[max_num_unique_values];
  this->counter_per_class = new size_t[num_classes * max_num_unique_values];
}

bool ProbabilitySplittingRule::findBestSplit(size_t nodeID,
                                             const std::vector<size_t>& possible_split_varIDs,
                                             std::unordered_map<size_t, double> responses_by_sampleID,
                                             std::vector<std::vector<size_t>> &sampleIDs,
                                             std::vector<size_t> &split_varIDs,
                                             std::vector<double> &split_values) {
  size_t num_samples_node = sampleIDs[nodeID].size();
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;

  size_t* class_counts = new size_t[num_classes]();
  // Compute overall class counts
  for (size_t i = 0; i < num_samples_node; ++i) {
    size_t sampleID = sampleIDs[nodeID][i];
    uint sample_classID = (uint) round(responses_by_sampleID[sampleID]);
    ++class_counts[sample_classID];
  }

  // For all possible split variables
  for (auto &varID : possible_split_varIDs) {
    // Use faster method for both cases
    double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
    if (q < Q_THRESHOLD) {
      findBestSplitValueSmallQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
                               best_decrease, responses_by_sampleID, sampleIDs);
    } else {
      findBestSplitValueLargeQ(nodeID, varID, num_classes, class_counts, num_samples_node, best_value, best_varID,
                               best_decrease, responses_by_sampleID, sampleIDs);
    }
  }

  delete[] class_counts;

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;
  return false;
}

void ProbabilitySplittingRule::findBestSplitValueSmallQ(size_t nodeID, size_t varID,size_t num_classes,
                                                        size_t *class_counts,
                                                        size_t num_samples_node, double &best_value, size_t &best_varID,
                                                        double &best_decrease,
                                                        std::unordered_map<size_t, double> responses_by_sampleID,
                                                        std::vector<std::vector<size_t>> &sampleIDs) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs[nodeID], varID);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Remove largest value because no split possible
  possible_split_values.pop_back();

  // Initialize with 0, if not in memory efficient mode, use pre-allocated space
  size_t num_splits = possible_split_values.size();
  size_t *class_counts_right = counter_per_class;
  size_t *n_right = counter;

  std::fill(class_counts_right, class_counts_right + num_splits * num_classes, 0);
  std::fill(n_right, n_right + num_splits, 0);

  // Count samples in right child per class and possbile split
  for (auto& sampleID : sampleIDs[nodeID]) {
    double value = data->get(sampleID, varID);
    uint sample_classID = responses_by_sampleID[sampleID];

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        ++class_counts_right[i * num_classes + sample_classID];
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

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      size_t class_count_right = class_counts_right[i * num_classes + j];
      size_t class_count_left = class_counts[j] - class_count_right;

      sum_right += class_count_right * class_count_right;
      sum_left += class_count_left * class_count_left;
    }

    // Decrease of impurity
    double decrease = sum_left / (double) n_left + sum_right / (double) n_right[i];

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = possible_split_values[i];
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void ProbabilitySplittingRule::findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes,
                                                        size_t *class_counts,
                                                        size_t num_samples_node, double &best_value, size_t &best_varID,
                                                        double &best_decrease,
                                                        std::unordered_map<size_t, double> responses_by_sampleID,
                                                        std::vector<std::vector<size_t>> &sampleIDs) {
  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill(counter_per_class, counter_per_class + num_unique * num_classes, 0);
  std::fill(counter, counter + num_unique, 0);

  // Count values
  for (auto &sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);
    size_t classID = responses_by_sampleID[sampleID];

    ++counter[index];
    ++counter_per_class[index * num_classes + classID];
  }

  size_t n_left = 0;
  size_t* class_counts_left = new size_t[num_classes]();

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {

    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    // Sum of squares
    double sum_left = 0;
    double sum_right = 0;
    for (size_t j = 0; j < num_classes; ++j) {
      class_counts_left[j] += counter_per_class[i * num_classes + j];
      size_t class_count_right = class_counts[j] - class_counts_left[j];

      sum_left += class_counts_left[j] * class_counts_left[j];
      sum_right += class_count_right * class_count_right;
    }

    // Decrease of impurity
    double decrease = sum_right / (double) n_right + sum_left / (double) n_left;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = data->getUniqueDataValue(varID, i);
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}
