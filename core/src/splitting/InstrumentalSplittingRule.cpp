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

#include "InstrumentalSplittingRule.h"

namespace grf {

InstrumentalSplittingRule::InstrumentalSplittingRule(size_t max_num_unique_values,
                                                     uint min_node_size,
                                                     double alpha,
                                                     double imbalance_penalty):
    min_node_size(min_node_size),
    alpha(alpha),
    imbalance_penalty(imbalance_penalty) {
  this->counter = new size_t[max_num_unique_values];
  this->sums = new double[max_num_unique_values];
  this->num_small_z = new size_t[max_num_unique_values];
  this->sums_z = new double[max_num_unique_values];
  this->sums_z_squared = new double[max_num_unique_values];
}

InstrumentalSplittingRule::~InstrumentalSplittingRule() {
  if (counter != nullptr) {
    delete[] counter;
  }
  if (sums != nullptr) {
    delete[] sums;
  }
  if (sums_z != nullptr) {
    delete[] sums_z;
  }
  if (sums_z_squared != nullptr) {
    delete[] sums_z_squared;
  }
  if (num_small_z != nullptr) {
    delete[] num_small_z;
  }
}

bool InstrumentalSplittingRule::find_best_split(const Data& data,
                                                size_t node,
                                                const std::vector<size_t>& possible_split_vars,
                                                const std::vector<double>& responses_by_sample,
                                                const std::vector<std::vector<size_t>>& samples,
                                                std::vector<size_t>& split_vars,
                                                std::vector<double>& split_values) {
  size_t num_samples = samples[node].size();

  // Precompute relevant quantities for this node.
  double sum_node = 0.0;
  double sum_node_z = 0.0;
  double sum_node_z_squared = 0.0;
  for (auto& sample : samples[node]) {
    sum_node += responses_by_sample[sample];

    double z = data.get_instrument(sample);
    sum_node_z += z;
    sum_node_z_squared += z * z;
  }

  double size_node = sum_node_z_squared - sum_node_z * sum_node_z / (double) num_samples;
  double min_child_size = size_node * alpha;

  double mean_z_node = sum_node_z / num_samples;
  size_t num_node_small_z = 0;
  for (auto& sample : samples[node]) {
    double z = data.get_instrument(sample);
    if (z < mean_z_node) {
      num_node_small_z++;
    }
  }

  // Initialize the variables to track the best split variable.
  size_t best_var = 0;
  double best_value = 0;
  double best_decrease = 0.0;

  for (auto& var : possible_split_vars) {
    // Use faster method for both cases
    double q = (double) num_samples / (double) data.get_num_unique_data_values(var);
    if (q < Q_THRESHOLD) {
      find_best_split_value_small_q(data, node, var, num_samples, sum_node, mean_z_node, num_node_small_z,
                                    sum_node_z, sum_node_z_squared, min_child_size, best_value,
                                    best_var, best_decrease, responses_by_sample, samples);
    } else {
      find_best_split_value_large_q(data, node, var, num_samples, sum_node, mean_z_node, num_node_small_z,
                                    sum_node_z, sum_node_z_squared, min_child_size, best_value,
                                    best_var, best_decrease, responses_by_sample, samples);
    }
  }

  // Stop if no good split found
  if (best_decrease <= 0.0) {
    return true;
  }

  // Save best values
  split_vars[node] = best_var;
  split_values[node] = best_value;
  return false;
}

void InstrumentalSplittingRule::find_best_split_value_small_q(const Data& data,
                                                              size_t node, size_t var,
                                                              size_t num_samples,
                                                              double sum_node,
                                                              double mean_node_z,
                                                              size_t num_node_small_z,
                                                              double sum_node_z,
                                                              double sum_node_z_squared,
                                                              double min_child_size,
                                                              double& best_value,
                                                              size_t& best_var,
                                                              double& best_decrease,
                                                              const std::vector<double>& responses_by_sample,
                                                              const std::vector<std::vector<size_t>>& samples) {
  std::vector<double> possible_split_values;
  data.get_all_values(possible_split_values, samples[node], var);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Remove largest value because no split possible
  possible_split_values.pop_back();

  // Initialize with 0m if not in memory efficient mode, use pre-allocated space
  size_t num_splits = possible_split_values.size();
  double* sums_right = sums;
  size_t* n_right = counter;

  std::fill(sums_right, sums_right + num_splits, 0);
  std::fill(n_right, n_right + num_splits, 0);
  std::fill(num_small_z, num_small_z + num_splits, 0);
  std::fill(sums_z, sums_z + num_splits, 0);
  std::fill(sums_z_squared, sums_z_squared + num_splits, 0);

  // Sum in right child and possible split
  for (auto& sample : samples[node]) {
    double value = data.get(sample, var);
    double label = responses_by_sample[sample];
    double z = data.get_instrument(sample);

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        sums_right[i] += label;
        sums_z[i] += z;
        sums_z_squared[i] += z * z;
        if (z < mean_node_z) {
          ++num_small_z[i];
        }
      } else {
        break;
      }
    }
  }

  // Compute decrease of impurity for each possible split.
  for (size_t i = 0; i < num_splits; ++i) {
    size_t n_left = num_samples - n_right[i];

    // Stop if the right child does not contain enough z
    // values below and above the parent's mean.
    size_t num_right_small_z = num_small_z[i];
    size_t num_right_large_z = n_right[i] - num_right_small_z;
    if (num_right_small_z < min_node_size || num_right_large_z < min_node_size) {
      break;
    }

    // Skip this split if the left child does not contain enough
    // z values below and above the parent's mean.
    size_t num_left_small_z = num_node_small_z - num_right_small_z;
    size_t num_left_large_z = n_left - num_left_small_z;
    if (num_left_small_z < min_node_size || num_left_large_z < min_node_size) {
      continue;
    }

    // Calculate relevant quantities for the right child.
    double sum_right = sums_right[i];
    double sum_right_z = sums_z[i];
    double sum_right_z_squared = sums_z_squared[i];
    double size_right = sum_right_z_squared - sum_right_z * sum_right_z / (double) n_right[i];

    // Skip this split if the right child's variance is too small.
    if (size_right < min_child_size || (imbalance_penalty > 0.0 && size_right == 0)) {
      continue;
    }

    // Calculate relevant quantities for the left child.
    double sum_left = sum_node - sum_right;
    double sum_left_z = sum_node_z - sum_right_z;
    double sum_left_z_squared = sum_node_z_squared - sum_right_z_squared;
    double size_left = sum_left_z_squared - sum_left_z * sum_left_z / (double) n_left;

    // Skip this split if the left child's variance is too small.
    if (size_left < min_child_size || (imbalance_penalty > 0.0 && size_left == 0)) {
      continue;
    }

    // Calculate the decrease in impurity.
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right[i];
    decrease -= imbalance_penalty * (1.0 / size_left + 1.0 / size_right);

    // Save this split if it is the best seen so fa.
    if (decrease > best_decrease) {
      best_value = possible_split_values[i];
      best_var = var;
      best_decrease = decrease;
    }
  }
}

void InstrumentalSplittingRule::find_best_split_value_large_q(const Data& data,
                                                              size_t node,
                                                              size_t var,
                                                              size_t num_samples,
                                                              double sum_node,
                                                              double mean_node_z,
                                                              size_t num_node_small_z,
                                                              double sum_node_z,
                                                              double sum_node_z_squared,
                                                              double min_child_size,
                                                              double& best_value,
                                                              size_t& best_var,
                                                              double& best_decrease,
                                                              const std::vector<double>& responses_by_sample,
                                                              const std::vector<std::vector<size_t>>& samples) {
  // Set counters to 0
  size_t num_unique = data.get_num_unique_data_values(var);
  std::fill(counter, counter + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);
  std::fill(num_small_z, num_small_z + num_unique, 0);
  std::fill(sums_z, sums_z + num_unique, 0);
  std::fill(sums_z_squared, sums_z_squared + num_unique, 0);

  for (auto& sample : samples[node]) {
    size_t i = data.get_index(sample, var);
    double z = data.get_instrument(sample);

    sums[i] += responses_by_sample[sample];
    ++counter[i];

    sums_z[i] += z;
    sums_z_squared[i] += z * z;
    if (z < mean_node_z) {
      ++num_small_z[i];
    }
  }

  size_t n_left = 0;
  double sum_left = 0;
  double sum_left_z = 0.0;
  double sum_left_z_squared = 0.0;
  size_t num_left_small_z = 0;

  // Compute decrease of impurity for each possible split.
  for (size_t i = 0; i < num_unique - 1; ++i) {
    n_left += counter[i];
    num_left_small_z += num_small_z[i];
    sum_left += sums[i];
    sum_left_z += sums_z[i];
    sum_left_z_squared += sums_z_squared[i];

    // Skip this split if the left child does not contain enough
    // z values below and above the parent's mean.
    size_t num_left_large_z = n_left - num_left_small_z;
    if (num_left_small_z < min_node_size || num_left_large_z < min_node_size) {
      continue;
    }

    // Stop if the right child does not contain enough z values below
    // and above the parent's mean.
    size_t n_right = num_samples - n_left;
    size_t num_right_small_z = num_node_small_z - num_left_small_z;
    size_t num_right_large_z = n_right - num_right_small_z;
    if (num_right_small_z < min_node_size || num_right_large_z < min_node_size) {
      break;
    }

    // Calculate relevant quantities for the left child.
    double size_left = sum_left_z_squared - sum_left_z * sum_left_z / (double) n_left;

    // Skip this split if the left child's variance is too small.
    if (size_left < min_child_size || (imbalance_penalty > 0.0 && size_left == 0)) {
      continue;
    }

    // Calculate relevant quantities for the left child.
    double sum_right_z_squared = sum_node_z_squared - sum_left_z_squared;
    double sum_right_z = sum_node_z - sum_left_z;
    double size_right = sum_right_z_squared - sum_right_z * sum_right_z / (double) n_right;

    // Skip this split if the right child's variance is too small.
    if (size_right < min_child_size || (imbalance_penalty > 0.0 && size_right == 0)) {
      continue;
    }

    // Calculate the decrease in impurity.
    double sum_right = sum_node - sum_left;
    double decrease = sum_left * sum_left / (double) n_left + sum_right * sum_right / (double) n_right;

    // Penalize splits that are too close to the edges of the data.
    decrease -= imbalance_penalty * (1.0 / size_left + 1.0 / size_right);

    // Save this split if it is the best seen so fa.
    if (decrease > best_decrease) {
      best_value = data.get_unique_data_value(var, i);
      best_var = var;
      best_decrease = decrease;
    }
  }
}

} // namespace grf
