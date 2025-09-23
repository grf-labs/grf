/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#include "AcceleratedSurvivalSplittingRule.h"

namespace grf {

AcceleratedSurvivalSplittingRule::AcceleratedSurvivalSplittingRule(double alpha):
    alpha(alpha) {
}

bool AcceleratedSurvivalSplittingRule::find_best_split(const Data& data,
                                            size_t node,
                                            const std::vector<size_t>& possible_split_vars,
                                            const Eigen::ArrayXXd& responses_by_sample,
                                            const std::vector<std::vector<size_t>>& samples_by_node,
                                            std::vector<size_t>& split_vars,
                                            std::vector<double>& split_values,
                                            std::vector<bool>& send_missing_left) {
  const std::vector<size_t>& samples = samples_by_node[node];

  // The splitting rule output
  double best_value = 0;
  size_t best_var = 0;
  bool best_send_missing_left = true;
  double best_logrank = 0;

  find_best_split_internal(data, possible_split_vars, responses_by_sample, samples,
                           best_value, best_var, best_send_missing_left, best_logrank);

  // Stop if no good split found
  if (best_logrank <= 0.0) {
    return true;
  }

  // Save best values
  split_vars[node] = best_var;
  split_values[node] = best_value;
  send_missing_left[node] = best_send_missing_left;
  return false;
}

void AcceleratedSurvivalSplittingRule::find_best_split_internal(const Data& data,
                                                     const std::vector<size_t>& possible_split_vars,
                                                     const Eigen::ArrayXXd& responses_by_sample,
                                                     const std::vector<size_t>& samples,
                                                     double& best_value,
                                                     size_t& best_var,
                                                     bool& best_send_missing_left,
                                                     double& best_logrank) {
  size_t size_node = samples.size();
  size_t min_child_size = std::max<size_t>(static_cast<size_t>(std::ceil(size_node * alpha)), 1uL);

  // Get the failure values t1, ..., tm in this node
  std::vector<double> failure_values;
  for (auto& sample : samples) {
    if (data.is_failure(sample)) {
      failure_values.push_back(responses_by_sample(sample, 0));
    }
  }

  size_t num_failures_node = failure_values.size();
  std::sort(failure_values.begin(), failure_values.end());
  failure_values.erase(std::unique(failure_values.begin(), failure_values.end()), failure_values.end());

  // The number of unique failure values in this node
  size_t num_failures = failure_values.size();
  // If there are no failures or only one failure time there is nothing to do.
  if (num_failures <= 1) {
    return;
  }

  // The number of failures at each time in the parent node. Entry 0 will be zero.
  // (Entry 0 is for time k < t1)
  std::vector<double> count_failure(num_failures + 1);
  // The number of censored observations at each time in the parent node.
  std::vector<double> count_censor(num_failures + 1);
  // The number of samples in the parent node at risk at each time point, i.e. the count of observations
  // with observed time greater than or equal to the given failure time. Entry 0 will be equal to the number
  // of samples (and the entries will always be monotonically decreasing)
  std::vector<double> at_risk(num_failures + 1);
  at_risk[0] = static_cast<double>(size_node);

  // allocating an N-sized (full data set size) array is faster than a hash table
  std::vector<size_t> relabeled_failures(data.get_num_rows());

  std::vector<double> numerator_weights(num_failures + 1);
  std::vector<double> cumsum_weights(num_failures + 1);

  // Relabel the failure values to range from 0 to the number of failures in this node
  for (auto& sample : samples) {
    double failure_value = responses_by_sample(sample, 0);
    size_t new_failure_value = std::upper_bound(failure_values.begin(), failure_values.end(),
                                                failure_value) - failure_values.begin();
    relabeled_failures[sample] = new_failure_value;
    if (data.is_failure(sample)) {
      ++count_failure[new_failure_value];
    } else {
      ++count_censor[new_failure_value];
    }
  }

  for (size_t time = 1; time < num_failures + 1; time++) {
    at_risk[time] = at_risk[time - 1] - count_failure[time - 1] - count_censor[time - 1];
    double Yk = at_risk[time];
    double dk = count_failure[time];
    numerator_weights[time] = dk / Yk;
  }

  for (size_t time = 1; time < num_failures + 1; time++) {
    cumsum_weights[time] = cumsum_weights[time - 1] + numerator_weights[time];
  }

  double gamma_node = 0;
  for (auto& sample : samples) {
    size_t sample_time = relabeled_failures[sample];
    gamma_node += cumsum_weights[sample_time];
  }

  for (auto& var : possible_split_vars) {
    find_best_split_value(data, var, size_node, min_child_size, num_failures_node, num_failures,
                          best_value, best_var, best_logrank, best_send_missing_left, samples,
                          relabeled_failures, cumsum_weights, gamma_node);
  }
}

void AcceleratedSurvivalSplittingRule::find_best_split_value(const Data& data,
                                                  size_t var,
                                                  size_t size_node,
                                                  size_t min_child_size,
                                                  size_t num_failures_node,
                                                  size_t num_failures,
                                                  double& best_value,
                                                  size_t& best_var,
                                                  double& best_logrank,
                                                  bool& best_send_missing_left,
                                                  const std::vector<size_t>& samples,
                                                  const std::vector<size_t>& relabeled_failures,
                                                  const std::vector<double>& cumsum_weights,
                                                  double gamma_node) {
  // possible_split_values contains all the unique split values for this variable in increasing order
  // sorted_samples contains the samples in this node in increasing order
  // if there are missing values, these are placed first
  // (if all Xij's are continuous, these two vectors have the same length)
  std::vector<double> possible_split_values;
  std::vector<size_t> sorted_samples;
  data.get_all_values(possible_split_values, sorted_samples, samples, var);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  size_t n_missing = 0;
  size_t num_failures_missing = 0;
  double numerator_missing = 0;
  double E1_missing = 0;

  // Loop through all samples to scan for missing values
  for (size_t i = 0; i < size_node - 1; i++) {
    size_t sample = sorted_samples[i];
    double sample_value = data.get(sample, var);
    if (std::isnan(sample_value)) {
      size_t sample_time = relabeled_failures[sample];
      double delta = data.is_failure(sample) ? 1.0 : 0.0;
      double gammai = cumsum_weights[sample_time];
      numerator_missing += delta - gammai;
      E1_missing += gammai;
      ++n_missing;
      if (delta > 0) {
        ++num_failures_missing;
      }
    }
  }

  size_t num_splits = possible_split_values.size() - 1;
  size_t num_failures_left = num_failures_missing;
  size_t split_index = 0;
  size_t start_sample = n_missing > 0 ? n_missing - 1 : 0;
  double numerator = numerator_missing;
  double E1 = E1_missing;
  double E2 = gamma_node - E1;

  for (bool send_left : {true, false}) {
    if (!send_left) {
     // A normal split with no NaNs, so we can stop early.
     if (n_missing == 0) {
       break;
     }
     // Else, send all missing right
     num_failures_left = 0;
     numerator = 0;
     E1 = 0;
     // Not necessary to evaluate splitting on NaN when sending right.
     split_index = 1;
     start_sample = n_missing;
    }

    for (size_t i = start_sample; i < size_node - 1; i++) {
      size_t sample = sorted_samples[i];
      size_t next_sample = sorted_samples[i + 1];
      double sample_value = data.get(sample, var);
      double next_sample_value = data.get(next_sample, var);
      size_t sample_time = relabeled_failures[sample];
      double delta = data.is_failure(sample) ? 1.0 : 0.0;

      // If there are missing values, we evaluate splitting on NaN when send_left is true
      // and i = n_missing - 1, which is why we need to check for missing below.
      bool split_on_missing = std::isnan(sample_value);

      if (!split_on_missing) {
        double gammai = cumsum_weights[sample_time];
        numerator += delta - gammai;
        E1 += gammai;
        E2 = gamma_node - E1;
        if (delta > 0) {
          ++num_failures_left;
        }
      }

      // Skip this split if one child is too small (i.e. too few failures)
      if (num_failures_left < min_child_size) {
        if (sample_value != next_sample_value) {
          ++split_index;
        }
        continue;
      }

      // Stop if the right child is too small.
      size_t num_failures_right = num_failures_node - num_failures_left;
      if (num_failures_right < min_child_size) {
        break;
      }

      // If the next sample value is different we can evaluate a split here
      if (sample_value != next_sample_value) {
        // The numerator is identical to the one in SurvivalSplittingRule, but calculated in O(1) updates
        // The denominator is an approximation, but which can be updated in O(1) (details are in arXiv TODO)
        double logrank = 0;
        if (E1 > 0 && E2 > 0) {
          double denominator = 1 / E1 + 1 / E2;
          logrank = numerator * numerator * denominator;
        }

        if (logrank > best_logrank) {
          best_value = possible_split_values[split_index];
          best_var = var;
          best_logrank = logrank;
          best_send_missing_left = send_left;
        }
        ++split_index;
      }

      if (split_index == num_splits) {
        break;
      }
    }
  }
}

} // namespace grf
