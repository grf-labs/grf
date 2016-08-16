
#include <set>
#include <unordered_map>
#include "TreeCausal.h"

TreeCausal::TreeCausal(size_t treatment_varID) : treatment_varID(treatment_varID),
                                                 udist(std::uniform_int_distribution<uint>()) {}

bool TreeCausal::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    split_values[nodeID] = getTerminalNodePrediction(sampleIDs[nodeID]);
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
    split_values[nodeID] = getTerminalNodePrediction(sampleIDs[nodeID]);
    return true;
  }

  std::unordered_map<size_t, double> responses_by_sampleIDs = relabelResponses(sampleIDs[nodeID]);
  bool stop = findBestSplit(nodeID, possible_split_varIDs, responses_by_sampleIDs);
  if (stop) {
    split_values[nodeID] = getTerminalNodePrediction(sampleIDs[nodeID]);
    return true;
  }

  return false;
}

std::unordered_map<size_t, double> TreeCausal::relabelResponses(std::vector<size_t>& node_sampleIDs) {
  std::unordered_map<size_t, double> new_responses_by_sampleID;

  size_t num_treated_responses = 0;
  for (size_t sampleID : node_sampleIDs) {
    num_treated_responses += data->get(sampleID, treatment_varID);
  }
  double average_treatment = (double) num_treated_responses / node_sampleIDs.size();

  std::unordered_map<bool, double> average_responses_by_treatment = calculateAverageResponses(node_sampleIDs);

  for (size_t sampleID : node_sampleIDs) {
    double treatment = data->get(sampleID, treatment_varID);
    double response = data->get(sampleID, dependent_varID);

    double new_response;
    if (equalDoubles(0, treatment)) {
      new_response = -(response - average_responses_by_treatment[false]) / (1 - average_treatment);
    } else if (equalDoubles(1, treatment)) {
      new_response = (response - average_responses_by_treatment[true]) / average_treatment;
    } else {
      throw std::runtime_error("Treatment values must be either 0 or 1.");
    }
    new_responses_by_sampleID[sampleID] = new_response;
  }
  return new_responses_by_sampleID;
}

double TreeCausal::getTerminalNodePrediction(std::vector<size_t>& node_sampleIDs) {
  std::unordered_map<bool, double> average_responses_by_treatment = calculateAverageResponses(node_sampleIDs);
  return average_responses_by_treatment[true] - average_responses_by_treatment[false];
}

std::unordered_map<bool, double> TreeCausal::calculateAverageResponses(std::vector<size_t>& node_sampleIDs) {
  std::unordered_map<bool, int> num_responses_by_treatment;
  std::unordered_map<bool, double> sum_responses_by_treatment;

  for (size_t sampleID : node_sampleIDs) {
    double treatment = data->get(sampleID, treatment_varID);
    double response = data->get(sampleID, dependent_varID);
    if (equalDoubles(0, treatment)) {
      sum_responses_by_treatment[false] += response;
      num_responses_by_treatment[false]++;
    } else if (equalDoubles(1, treatment)) {
      sum_responses_by_treatment[true] += response;
      num_responses_by_treatment[true]++;
    } else {
      throw std::runtime_error("Treatment values must be either 0 or 1.");
    }
  }

  if (num_responses_by_treatment[true] <= 0) {
    throw std::runtime_error("Encountered a pure untreated node.");
  }

  if (num_responses_by_treatment[false] <= 0) {
    throw std::runtime_error("Encountered a pure treated node.");
  }

  std::unordered_map<bool, double> average_responses_by_treatment;
  for (auto it = sum_responses_by_treatment.begin(); it != sum_responses_by_treatment.end(); ++it) {
    bool treatment = it->first;
    average_responses_by_treatment[treatment] = it->second / num_responses_by_treatment[treatment];
  }
  return average_responses_by_treatment;
};

bool TreeCausal::equalDoubles(double first, double second) {
 return std::abs(first - second) < std::numeric_limits<double>::epsilon();
}

// NOTE: This rest of this file is largely copied and pasted from TreeRegression, with small
// modifications to avoid a s plit that would result in homogeneous treatment values.

void TreeCausal::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                          double &best_value, size_t &best_varID, double &best_decrease,
                                          std::unordered_map<size_t, double> &responses_by_sampleID) {

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
  std::vector<double> treatment_sums_right(num_splits);
  double total_treatment = 0.0;

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
    double response = responses_by_sampleID[sampleID];

    double treatment = data->get(sampleID, treatment_varID);
    total_treatment += treatment;

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        sums_right[i] += response;
        treatment_sums_right[i] += treatment;
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

    double treatment_sum_right = treatment_sums_right[i];
    if (equalDoubles(treatment_sum_right, 0)
        || equalDoubles(treatment_sum_right, n_right[i])
        || equalDoubles(treatment_sum_right, total_treatment)
        || equalDoubles(treatment_sum_right + n_left, total_treatment)) {
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

void TreeCausal::findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                          double &best_value, size_t &best_varID, double &best_decrease,
                                          std::unordered_map<size_t, double> &responses_by_sampleID) {

  // Set counters to 0
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill(counter, counter + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);

  std::vector<double> treatment_sums(num_unique);
  double total_treatment = 0.0;

  for (auto& sampleID : sampleIDs[nodeID]) {
    size_t index = data->getIndex(sampleID, varID);

    sums[index] += responses_by_sampleID[sampleID];
    ++counter[index];
    treatment_sums[index] += data->get(sampleID, treatment_varID);
    total_treatment += data->get(sampleID, treatment_varID);
  }

  size_t n_left = 0;
  double sum_left = 0;
  double treatment_left = 0.0;

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {

    // Stop if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];
    sum_left += sums[i];
    treatment_left += treatment_sums[i];

    // Stop if right child empty
    size_t n_right = num_samples_node - n_left;
    if (n_right == 0) {
      break;
    }

    if (equalDoubles(treatment_left, 0)
        || equalDoubles(treatment_left, n_left)
        || equalDoubles(treatment_left, total_treatment)
        || equalDoubles(treatment_left + n_right, total_treatment)) {
      continue;
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