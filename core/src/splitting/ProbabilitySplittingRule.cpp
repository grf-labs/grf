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
#include <cmath>
#include <unordered_map>

#include "ProbabilitySplittingRule.h"

ProbabilitySplittingRule::ProbabilitySplittingRule(Data* data, double alpha, size_t num_classes) {
  this->data = data;
  this->alpha = alpha;
  this->num_classes = num_classes;

  size_t max_num_unique_values = data->get_max_num_unique_values();
  this->counter = new size_t[max_num_unique_values];
  this->counter_per_class = new size_t[num_classes * max_num_unique_values];
}

ProbabilitySplittingRule::~ProbabilitySplittingRule() {
  if (counter != 0) {
    delete[] counter;
  }
  if (counter_per_class != 0) {
    delete[] counter_per_class;
  }
}

bool ProbabilitySplittingRule::find_best_split(size_t node,
                                               const std::vector<size_t>& possible_split_vars,
                                               const std::unordered_map<size_t, double>& responses_by_sample,
                                               const std::vector<std::vector<size_t>>& samples,
                                               std::vector<size_t>& split_vars,
                                               std::vector<double>& split_values) {
  size_t num_samples_node = samples[node].size();
  size_t min_child_samples = std::max<size_t>(std::ceil(num_samples_node * alpha), 1uL);

  double best_decrease = -1;
  size_t best_var = 0;
  double best_value = 0;

  size_t* class_counts = new size_t[num_classes]();
  // Compute overall class counts
  for (size_t i = 0; i < num_samples_node; ++i) {
    size_t sample = samples[node][i];
    uint sample_class = (uint) std::round(responses_by_sample.at(sample));
    ++class_counts[sample_class];
  }

  // For all possible split variables
  for (size_t var : possible_split_vars) {
    // Use faster method for both cases
    double q = (double) num_samples_node / (double) data->get_num_unique_data_values(var);
    if (q < Q_THRESHOLD) {
      find_best_split_value_small_q(node, var, num_classes, class_counts, num_samples_node, min_child_samples,
                                    best_value, best_var, best_decrease, responses_by_sample, samples);
    } else {
      find_best_split_value_large_q(node, var, num_classes, class_counts, num_samples_node, min_child_samples,
                                    best_value, best_var, best_decrease, responses_by_sample, samples);
    }
  }

  delete[] class_counts;

  // Stop if no good split found
  if (best_decrease < 0) {
    return true;
  }

  // Save best values
  split_vars[node] = best_var;
  split_values[node] = best_value;
  return false;
}

void ProbabilitySplittingRule::find_best_split_value_small_q(size_t node, size_t var, size_t num_classes,
                                                             size_t* class_counts,
                                                             size_t num_samples_node,
                                                             size_t min_child_samples,
                                                             double& best_value,
                                                             size_t& best_var,
                                                             double& best_decrease,
                                                             const std::unordered_map<size_t, double>& labels_by_sample,
                                                             const std::vector<std::vector<size_t>>& samples) {

  // Create possible split values
  std::vector<double> possible_split_values;
  data->get_all_values(possible_split_values, samples[node], var);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  // Remove largest value because no split possible
  possible_split_values.pop_back();

  // Initialize with 0, if not in memory efficient mode, use pre-allocated space
  size_t num_splits = possible_split_values.size();
  size_t* class_counts_right = counter_per_class;
  size_t* n_right = counter;

  std::fill(class_counts_right, class_counts_right + num_splits * num_classes, 0);
  std::fill(n_right, n_right + num_splits, 0);

  // Count samples in right child per class and possbile split
  for (auto& sample : samples[node]) {
    double value = data->get(sample, var);
    uint sample_class = labels_by_sample.at(sample);

    // Count samples until split_value reached
    for (size_t i = 0; i < num_splits; ++i) {
      if (value > possible_split_values[i]) {
        ++n_right[i];
        ++class_counts_right[i * num_classes + sample_class];
      } else {
        break;
      }
    }
  }

  // Compute decrease of impurity for each possible split
  for (size_t i = 0; i < num_splits; ++i) {

    // Skip this split if one child is too small.
    size_t n_left = num_samples_node - n_right[i];
    if (n_left < min_child_samples || n_right[i] < min_child_samples) {
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
      best_var = var;
      best_decrease = decrease;
    }
  }
}

void ProbabilitySplittingRule::find_best_split_value_large_q(size_t node, size_t var, size_t num_classes,
                                                             size_t* class_counts,
                                                             size_t num_samples_node,
                                                             size_t min_child_samples,
                                                             double& best_value,
                                                             size_t& best_var,
                                                             double& best_decrease,
                                                             const std::unordered_map<size_t, double>& responses_by_sample,
                                                             const std::vector<std::vector<size_t>>& samples) {
  // Set counters to 0
  size_t num_unique = data->get_num_unique_data_values(var);
  std::fill(counter_per_class, counter_per_class + num_unique * num_classes, 0);
  std::fill(counter, counter + num_unique, 0);

  // Count values
  for (auto& sample : samples[node]) {
    size_t index = data->get_index(sample, var);
    size_t sample_class = responses_by_sample.at(sample);

    ++counter[index];
    ++counter_per_class[index * num_classes + sample_class];
  }

  size_t n_left = 0;
  size_t* class_counts_left = new size_t[num_classes]();

  // Compute decrease of impurity for each split
  for (size_t i = 0; i < num_unique - 1; ++i) {
    // Skip this split if the left child is empty.
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i];

    // Stop if the right child is too small.
    size_t n_right = num_samples_node - n_left;
    if (n_right < min_child_samples) {
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

    // Skip to the next value if the left child is too small.
    if (n_left < min_child_samples) {
        continue;
    }

    // Decrease of impurity
    double decrease = sum_right / (double) n_right + sum_left / (double) n_left;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = data->get_unique_data_value(var, i);
      best_var = var;
      best_decrease = decrease;
    }
  }
  delete[] class_counts_left;
}
