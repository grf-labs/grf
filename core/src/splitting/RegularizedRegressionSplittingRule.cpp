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
#include <math.h>

#include "RegularizedRegressionSplittingRule.h"

RegularizedRegressionSplittingRule::RegularizedRegressionSplittingRule(Data* data,
                                                                       double lambda,
                                                                       bool downweight_penalty) {
  this->data = data;
  this->lambda = lambda;
  this->downweight_penalty = downweight_penalty;

  size_t max_num_unique_values = data->get_max_num_unique_values();
  this->counter = new size_t[max_num_unique_values];
  this->sums = new double[max_num_unique_values];
}

RegularizedRegressionSplittingRule::~RegularizedRegressionSplittingRule() {
  if (counter != 0) {
    delete[] counter;
  }
  if (sums != 0) {
    delete[] sums;
  }
}

bool RegularizedRegressionSplittingRule::find_best_split(size_t node,
                                                         const std::vector<size_t>& possible_split_vars,
                                                         const std::unordered_map<size_t, double>& labels_by_sample,
                                                         const std::vector<std::vector<size_t>>& samples,
                                                         std::vector<size_t>& split_vars,
                                                         std::vector<double>& split_values) {

  size_t num_samples_node = samples[node].size();
  size_t best_var = 0;
  double best_value = 0;

  // Compute sum of outcomes in this node.
  double sum_node = 0.0;
  for (auto& sample : samples[node]) {
    sum_node += labels_by_sample.at(sample);
  }

  // Compute the outcome variance within this node.
  double node_mean = sum_node / num_samples_node;
  double node_impurity = 0.0;
  for (auto& sample : samples[node]) {
    double difference = labels_by_sample.at(sample) - node_mean;
    node_impurity += difference * difference;
  }

  double initial_decrease = -2 * lambda * node_impurity;
  double best_decrease = initial_decrease;

  // For all possible split variables
  for (auto& var : possible_split_vars) {
    // Use faster method for both cases
    double q = (double) num_samples_node / (double) data->get_num_unique_data_values(var);
    if (q < Q_THRESHOLD) {
      find_best_split_value_small_q(node, var, sum_node, num_samples_node, node_impurity,
                                    best_value, best_var, best_decrease, labels_by_sample, samples);
    } else {
      find_best_split_value_large_q(node, var, sum_node, num_samples_node, node_impurity,
                                    best_value, best_var, best_decrease, labels_by_sample, samples);
    }
  }

  // Stop if no good split found
  if (best_decrease <= initial_decrease) {
    return true;
  }

  // Save best values
  split_vars[node] = best_var;
  split_values[node] = best_value;
  return false;
}

void RegularizedRegressionSplittingRule::find_best_split_value_small_q(size_t node, size_t var, double sum_node,
                                                                       size_t num_samples_node,
                                                                       double node_impurity,
                                                                       double& best_value, size_t& best_var,
                                                                       double& best_decrease,
                                                                       const std::unordered_map<size_t, double>& responses_by_sample,
                                                                       const std::vector<std::vector<size_t>>& samples) {
  std::vector<double> possible_split_values;
  data->get_all_values(possible_split_values, samples.at(node), var);

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

  // Sum in right child and possible split
  for (auto& sample : samples[node]) {
    double value = data->get(sample, var);
    double response = responses_by_sample.at(sample);

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

    // Penalize splits that are too close to the edges of the data.
    double penalty = lambda * (1.0 / n_left + 1.0 / n_right[i]) * node_impurity;
    if (downweight_penalty) {
      penalty *= sqrt(num_samples_node);
    }
    decrease -= penalty;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = possible_split_values[i];
      best_var = var;
      best_decrease = decrease;
    }
  }
}

void RegularizedRegressionSplittingRule::find_best_split_value_large_q(size_t node, size_t var, double sum_node,
                                                                       size_t num_samples_node,
                                                                       double node_impurity,
                                                                       double& best_value,
                                                                       size_t& best_var,
                                                                       double& best_decrease,
                                                                       const std::unordered_map<size_t, double>& responses_by_sample,
                                                                       const std::vector<std::vector<size_t>>& samples) {

  // Set counters to 0
  size_t num_unique = data->get_num_unique_data_values(var);
  std::fill(counter, counter + num_unique, 0);
  std::fill(sums, sums + num_unique, 0);

  for (auto& sample : samples[node]) {
    size_t index = data->get_index(sample, var);

    sums[index] += responses_by_sample.at(sample);
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

    // Penalize splits that are too close to the edges of the data.
    double penalty = lambda * (1.0 / n_left + 1.0 / n_right) * node_impurity;
    if (downweight_penalty) {
      penalty *= sqrt(num_samples_node);
    }
    decrease -= penalty;

    // If better than before, use this
    if (decrease > best_decrease) {
      best_value = data->get_unique_data_value(var, i);
      best_var = var;
      best_decrease = decrease;
    }
  }
}
