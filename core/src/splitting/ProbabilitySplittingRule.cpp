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

#include "ProbabilitySplittingRule.h"

namespace grf {

ProbabilitySplittingRule::ProbabilitySplittingRule(size_t max_num_unique_values,
                                                   size_t num_classes,
                                                   double alpha,
                                                   double imbalance_penalty) {
  this->num_classes = num_classes;

  this->alpha = alpha;
  this->imbalance_penalty = imbalance_penalty;

  this->counter = new size_t[max_num_unique_values];
  this->counter_per_class = new double[num_classes * max_num_unique_values];
}

ProbabilitySplittingRule::~ProbabilitySplittingRule() {
  if (counter != nullptr) {
    delete[] counter;
  }
  if (counter_per_class != nullptr) {
    delete[] counter_per_class;
  }
}

bool ProbabilitySplittingRule::find_best_split(const Data& data,
                                               size_t node,
                                               const std::vector<size_t>& possible_split_vars,
                                               const Eigen::ArrayXXd& responses_by_sample,
                                               const std::vector<std::vector<size_t>>& samples,
                                               std::vector<size_t>& split_vars,
                                               std::vector<double>& split_values,
                                               std::vector<bool>& send_missing_left) {
  size_t size_node = samples[node].size();
  size_t min_child_size = std::max<size_t>(static_cast<size_t>(std::ceil(size_node * alpha)), 1uL);

  double* class_counts = new double[num_classes]();
  for (size_t i = 0; i < size_node; ++i) {
    size_t sample = samples[node][i];
    uint sample_class = (uint) std::round(responses_by_sample(sample, 0));
    double sample_weight = data.get_weight(sample);
    class_counts[sample_class] += sample_weight;
  }

  // Initialize the variables to track the best split variable.
  size_t best_var = 0;
  double best_value = 0;
  double best_decrease = 0.0;
  bool best_send_missing_left = true;

  // For all possible split variables
  for (size_t var : possible_split_vars) {
    find_best_split_value(data, node, var, num_classes, class_counts, size_node, min_child_size,
                          best_value, best_var, best_decrease, best_send_missing_left, responses_by_sample, samples);
  }

  delete[] class_counts;

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

void ProbabilitySplittingRule::find_best_split_value(const Data& data,
                                                     size_t node, size_t var,
                                                     size_t num_classes,
                                                     double* class_counts,
                                                     size_t size_node,
                                                     size_t min_child_size,
                                                     double& best_value,
                                                     size_t& best_var,
                                                     double& best_decrease,
                                                     bool& best_send_missing_left,
                                                     const Eigen::ArrayXXd& responses_by_sample,
                                                     const std::vector<std::vector<size_t>>& samples) {
  std::vector<double> possible_split_values;
  std::vector<size_t> sorted_samples;
  data.get_all_values(possible_split_values, sorted_samples, samples[node], var);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  size_t num_splits = possible_split_values.size() - 1;

  std::fill(counter_per_class, counter_per_class + num_splits * num_classes, 0);
  std::fill(counter, counter + num_splits, 0);
  size_t n_missing = 0;
  double* class_counts_missing = new double[num_classes]();

  size_t split_index = 0;
  for (size_t i = 0; i < size_node - 1; i++) {
    size_t sample = sorted_samples[i];
    size_t next_sample = sorted_samples[i + 1];
    double sample_value = data.get(sample, var);
    uint sample_class = static_cast<uint>(responses_by_sample(sample, 0));
    double sample_weight = data.get_weight(sample);

    if (std::isnan(sample_value)) {
      class_counts_missing[sample_class] += sample_weight;
      ++n_missing;
    } else {
      ++counter[split_index];
      counter_per_class[split_index * num_classes + sample_class] += sample_weight;
    }

    double next_sample_value = data.get(next_sample, var);
    // if the next sample value is different, including the transition (..., NaN, Xij, ...)
    // then move on to the next bucket (all logical operators with NaN evaluates to false by default)
    if (sample_value != next_sample_value && !std::isnan(next_sample_value)) {
      ++split_index;
    }
  }

  size_t n_left = n_missing;
  double* class_counts_left = class_counts_missing;

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
      for (size_t cls = 0; cls < num_classes; ++cls) {
        class_counts_left[cls] = 0;
      }
    }

    for (size_t i = 0; i < num_splits; ++i) {
      // not necessary to evaluate sending right when splitting on NaN.
      if (i == 0 && !send_left) {
        continue;
      }

      n_left += counter[i];

      // Stop if the right child is too small.
      size_t n_right = size_node - n_left;
      if (n_right < min_child_size) {
        break;
      }

      // Sum of squares
      double sum_left = 0;
      double sum_right = 0;
      for (size_t cls = 0; cls < num_classes; ++cls) {
        class_counts_left[cls] += counter_per_class[i * num_classes + cls];
        double class_count_right = class_counts[cls] - class_counts_left[cls];

        sum_left += class_counts_left[cls] * class_counts_left[cls];
        sum_right += class_count_right * class_count_right;
      }

      // Skip to the next value if the left child is too small.
      if (n_left < min_child_size) {
          continue;
      }

      // Decrease of impurity
      double decrease = sum_right / (double) n_right + sum_left / (double) n_left;

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
  delete[] class_counts_missing;
}

} // namespace grf
