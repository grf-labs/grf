/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <algorithm>

#include "RegressionSplittingRule.h"

namespace grf {

RegressionSplittingRule::RegressionSplittingRule(size_t max_num_unique_values,
                                                 double alpha,
                                                 double imbalance_penalty):
    alpha(alpha),
    imbalance_penalty(imbalance_penalty) {
  this->counter = new size_t[max_num_unique_values];
  this->sums = new double[max_num_unique_values];
  this->weight_sums = new double[max_num_unique_values];
}

RegressionSplittingRule::~RegressionSplittingRule() {
  if (counter != nullptr) {
    delete[] counter;
  }
  if (sums != nullptr) {
    delete[] sums;
  }
  if (weight_sums != nullptr) {
    delete[] weight_sums;
  }
}

bool RegressionSplittingRule::find_best_split(const Data& data,
                                              size_t node,
                                              const std::vector<size_t>& possible_split_vars,
                                              const std::vector<double>& responses_by_sample,
                                              const std::vector<std::vector<size_t>>& samples,
                                              std::vector<size_t>& split_vars,
                                              std::vector<double>& split_values,
                                              std::vector<bool>& send_missing_left) {

  size_t size_node = samples[node].size();
  size_t min_child_size = std::max<size_t>(std::ceil(size_node * alpha), 1uL);

  // Precompute the sum of outcomes in this node.
  double sum_node = 0.0;
  double weight_sum_node = 0.0;
  for (auto& sample : samples[node]) {
    double sample_weight = data.get_weight(sample);
    weight_sum_node += sample_weight;
    sum_node += sample_weight * responses_by_sample[sample];
  }

  // Initialize the variables to track the best split variable.
  size_t best_var = 0;
  double best_value = 0;
  double best_decrease = 0.0;
  bool best_send_missing_left = true;

  // For all possible split variables
  for (auto& var : possible_split_vars) {
    // Use faster method for both cases
    double q = (double) size_node / (double) data.get_num_unique_data_values(var);
    if (q < Q_THRESHOLD || data.contains_nan()) {
      find_best_split_value_small_q(data, node, var, weight_sum_node, sum_node, size_node, min_child_size,
                                    best_value, best_var, best_decrease, best_send_missing_left, responses_by_sample, samples);
    } else {
      find_best_split_value_large_q(data, node, var, weight_sum_node, sum_node, size_node, min_child_size,
                                    best_value, best_var, best_decrease, responses_by_sample, samples);
    }
  }

  // Stop if no good split found
  if (best_decrease <= 0.0) {
    return true;
  }

  // Save best values
  split_vars[node] = best_var;
  split_values[node] = best_value;
  send_missing_left[node] = best_send_missing_left;
  return false;
}

void RegressionSplittingRule::find_best_split_value_small_q(const Data& data,
                                                            size_t node, size_t var,
                                                            double weight_sum_node,
                                                            double sum_node,
                                                            size_t size_node,
                                                            size_t min_child_size,
                                                            double& best_value, size_t& best_var,
                                                            double& best_decrease, bool& best_send_missing_left,
                                                            const std::vector<double>& responses_by_sample,
                                                            const std::vector<std::vector<size_t>>& samples) {
  // possible_split_values: the sorted unique split values. Length: num_splits (equal to size_node - 1 if all unique)
  // sorted_samples: the node samples in increasing order (may contain duplicated Xij). Length: size_node
  std::vector<double> possible_split_values;
  std::vector<size_t> sorted_samples;
  data.get_all_values(possible_split_values, sorted_samples, samples[node], var);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Initialize with 0m if not in memory efficient mode, use pre-allocated space
  size_t num_splits = possible_split_values.size() - 1; // -1: we do not split at the last value
  std::fill(weight_sums, weight_sums + num_splits, 0);
  std::fill(counter, counter + num_splits, 0);
  std::fill(sums, sums + num_splits, 0);
  size_t n_missing = 0;
  double weight_sum_missing = 0;
  double sum_missing = 0;

  // Fill counter and sums buckets
  size_t split_index = 0;
  for (size_t i = 0; i < size_node - 1; i++) {
    size_t sample = sorted_samples[i];
    size_t next_sample = sorted_samples[i + 1];
    double sample_value = data.get(sample, var);
    double response = responses_by_sample[sample];
    double sample_weight = data.get_weight(sample);

    if (std::isnan(sample_value)) {
      weight_sum_missing += sample_weight;
      sum_missing += sample_weight * response;
      ++n_missing;
    } else {
      weight_sums[split_index] += sample_weight;
      sums[split_index] += sample_weight * response;
      ++counter[split_index];
    }

    double next_sample_value = data.get(next_sample, var);
    // if the next sample value is different, including the transition (..., NaN, Xij, ...)
    // then move on to the next bucket (all logical operators with NaN evaluates to false by default)
    if (sample_value != next_sample_value && !std::isnan(next_sample_value)) {
      ++split_index;
    }
  }

  size_t n_left = n_missing;
  double weight_sum_left = weight_sum_missing;
  double sum_left = sum_missing;

  // Compute decrease of impurity for each possible split
  for (bool send_left : {true, false}) {
    if (!send_left) {
      // A normal split with no NaNs, so we can stop early.
      if (n_missing == 0) {
        break;
      }
      // It is not necessary to adjust n_right or sum_right as the the missing
      // part is included in the total sum.
      n_left = 0;
      weight_sum_left = 0;
      sum_left = 0;
    }

    for (size_t i = 0; i < num_splits; ++i) {
      // not necessary to evaluate sending right when splitting on NaN.
      if (i == 0 && !send_left) {
        continue;
      }

      n_left += counter[i];
      weight_sum_left += weight_sums[i];
      sum_left += sums[i];

      // Skip this split if one child is too small.
      if (n_left < min_child_size) {
        continue;
      }

      // Stop if the right child is too small.
      size_t n_right = size_node - n_left;
      if (n_right < min_child_size) {
        break;
      }

      double weight_sum_right = weight_sum_node - weight_sum_left;
      double sum_right = sum_node - sum_left;
      double decrease = sum_left * sum_left / weight_sum_left + sum_right * sum_right / weight_sum_right;

      // Penalize splits that are too close to the edges of the data.
      double penalty = imbalance_penalty * (1.0 / n_left + 1.0 / n_right);
      decrease -= penalty;


      // If better than before, use this
      if (decrease > best_decrease) {
        best_value = possible_split_values[i];
        best_var = var;
        best_decrease = decrease;
        best_send_missing_left = send_left;
      }
    }
  }
}

void RegressionSplittingRule::find_best_split_value_large_q(const Data& data,
                                                            size_t node,
                                                            size_t var,
                                                            double weight_sum_node,
                                                            double sum_node,
                                                            size_t size_node,
                                                            size_t min_child_size,
                                                            double& best_value,
                                                            size_t& best_var,
                                                            double& best_decrease,
                                                            const std::vector<double>& responses_by_sample,
                                                            const std::vector<std::vector<size_t>>& samples) {
  // Set counters to 0
  size_t num_unique = data.get_num_unique_data_values(var);
  std::fill(counter, counter + num_unique, 0);
  std::fill(weight_sums, weight_sums + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);
  for (auto& sample : samples[node]) {
    double sample_weight = data.get_weight(sample);
    size_t index = data.get_index(sample, var);
    weight_sums[index] += sample_weight;
    sums[index] += sample_weight * responses_by_sample[sample];
    ++counter[index];
  }

  size_t n_left = 0;
  double weight_sum_left = 0;
  double sum_left = 0;

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {
    // Continue if nothing here
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];
    weight_sum_left += weight_sums[i];
    sum_left += sums[i];

    // Skip to the next value if the left child is too small.
    if (n_left < min_child_size) {
      continue;
    }

    // Stop if the right child is too small.
    size_t n_right = size_node - n_left;
    if (n_right < min_child_size) {
      break;
    }

    double weight_sum_right = weight_sum_node - weight_sum_left;
    double sum_right = sum_node - sum_left;
    double decrease = sum_left * sum_left / weight_sum_left + sum_right * sum_right / weight_sum_right;

    // Penalize splits that are too close to the edges of the data.
    double penalty = imbalance_penalty * (1.0 / n_left + 1.0 / n_right);
    decrease -= penalty;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = data.get_unique_data_value(var, i);
      best_var = var;
      best_decrease = decrease;
    }
  }
}

} // namespace grf
