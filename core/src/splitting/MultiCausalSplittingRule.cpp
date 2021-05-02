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

#include "MultiCausalSplittingRule.h"

namespace grf {

MultiCausalSplittingRule::MultiCausalSplittingRule(size_t max_num_unique_values,
                                                   uint min_node_size,
                                                   double alpha,
                                                   double imbalance_penalty,
                                                   size_t response_length,
                                                   size_t num_treatments):
    min_node_size(min_node_size),
    alpha(alpha),
    imbalance_penalty(imbalance_penalty),
    response_length(response_length),
    num_treatments(num_treatments) {
  this->counter = new size_t[max_num_unique_values];
  this->weight_sums = new double[max_num_unique_values];
  this->sums = Eigen::ArrayXXd(max_num_unique_values, response_length);
  this->num_small_w = Eigen::ArrayXXi(max_num_unique_values, num_treatments);
  this->sums_w = Eigen::ArrayXXd(max_num_unique_values, num_treatments);
  this->sums_w_squared = Eigen::ArrayXXd(max_num_unique_values, num_treatments);
}

MultiCausalSplittingRule::~MultiCausalSplittingRule() {
  if (counter != nullptr) {
    delete[] counter;
  }
  if (weight_sums != nullptr) {
    delete[] weight_sums;
  }
}

bool MultiCausalSplittingRule::find_best_split(const Data& data,
                                               size_t node,
                                               const std::vector<size_t>& possible_split_vars,
                                               const Eigen::ArrayXXd& responses_by_sample,
                                               const std::vector<std::vector<size_t>>& samples,
                                               std::vector<size_t>& split_vars,
                                               std::vector<double>& split_values,
                                               std::vector<bool>& send_missing_left) {
  size_t num_samples = samples[node].size();

  // Precompute the sum of outcomes in this node.
  double weight_sum_node = 0.0;
  Eigen::ArrayXd sum_node = Eigen::ArrayXd::Zero(response_length);
  Eigen::ArrayXd sum_node_w = Eigen::ArrayXd::Zero(num_treatments);
  Eigen::ArrayXd sum_node_w_squared = Eigen::ArrayXd::Zero(num_treatments);
  // Allocate W-array and re-use to avoid expensive copy-inducing calls to `data.get_treatments`
  Eigen::ArrayXXd treatments = Eigen::ArrayXXd(num_samples, num_treatments);
  for (size_t i = 0; i < num_samples; i++) {
    size_t sample = samples[node][i];
    double sample_weight = data.get_weight(sample);
    weight_sum_node += sample_weight;
    sum_node += sample_weight * responses_by_sample.row(sample);
    treatments.row(i) = data.get_treatments(sample);

    sum_node_w += sample_weight * treatments.row(i);
    sum_node_w_squared += sample_weight * treatments.row(i).square();
  }

  Eigen::ArrayXd size_node = sum_node_w_squared - sum_node_w.square() / weight_sum_node;
  Eigen::ArrayXd min_child_size = size_node * alpha;

  Eigen::ArrayXd mean_w_node = sum_node_w / weight_sum_node;
  Eigen::ArrayXi num_node_small_w = Eigen::ArrayXi::Zero(num_treatments);
  for (size_t i = 0; i < num_samples; i++) {
    num_node_small_w += (treatments.row(i).transpose() < mean_w_node).cast<int>();
  }

  // Initialize the variables to track the best split variable.
  size_t best_var = 0;
  double best_value = 0;
  double best_decrease = 0.0;
  bool best_send_missing_left = true;

  // For all possible split variables
  for (auto& var : possible_split_vars) {
    find_best_split_value(data, node, var, num_samples, weight_sum_node, sum_node, mean_w_node, num_node_small_w,
                          sum_node_w, sum_node_w_squared, min_child_size, treatments, best_value,
                          best_var, best_decrease, best_send_missing_left, responses_by_sample, samples);
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

void MultiCausalSplittingRule::find_best_split_value(const Data& data,
                                                     size_t node,
                                                     size_t var,
                                                     size_t num_samples,
                                                     double weight_sum_node,
                                                     const Eigen::ArrayXd& sum_node,
                                                     const Eigen::ArrayXd& mean_node_w,
                                                     const Eigen::ArrayXi& num_node_small_w,
                                                     const Eigen::ArrayXd& sum_node_w,
                                                     const Eigen::ArrayXd& sum_node_w_squared,
                                                     const Eigen::ArrayXd& min_child_size,
                                                     const Eigen::ArrayXXd& treatments,
                                                     double& best_value,
                                                     size_t& best_var,
                                                     double& best_decrease,
                                                     bool& best_send_missing_left,
                                                     const Eigen::ArrayXXd& responses_by_sample,
                                                     const std::vector<std::vector<size_t>>& samples) {
  std::vector<double> possible_split_values;
  std::vector<size_t> sorted_samples;
  std::vector<size_t> index = data.get_all_values(possible_split_values, sorted_samples, samples[node], var);

  // Try next variable if all equal for this
  if (possible_split_values.size() < 2) {
    return;
  }

  size_t num_splits = possible_split_values.size() - 1;

  std::fill(counter, counter + num_splits, 0);
  std::fill(weight_sums, weight_sums + num_splits, 0);
  sums.topRows(num_splits).setZero();
  num_small_w.topRows(num_splits).setZero();
  sums_w.topRows(num_splits).setZero();
  sums_w_squared.topRows(num_splits).setZero();
  size_t n_missing = 0;
  double weight_sum_missing = 0;
  Eigen::ArrayXd sum_missing = Eigen::ArrayXd::Zero(response_length);
  Eigen::ArrayXd sum_w_missing = Eigen::ArrayXd::Zero(num_treatments);
  Eigen::ArrayXd sum_w_squared_missing = Eigen::ArrayXd::Zero(num_treatments);
  Eigen::ArrayXi num_small_w_missing = Eigen::ArrayXi::Zero(num_treatments);

  size_t split_index = 0;
  for (size_t i = 0; i < num_samples - 1; i++) {
    size_t sample = sorted_samples[i];
    size_t next_sample = sorted_samples[i + 1];
    size_t sort_index = index[i];
    double sample_value = data.get(sample, var);
    double sample_weight = data.get_weight(sample);

    if (std::isnan(sample_value)) {
      weight_sum_missing += sample_weight;
      sum_missing += sample_weight * responses_by_sample.row(sample);
      ++n_missing;

      sum_w_missing += sample_weight * treatments.row(sort_index);
      sum_w_squared_missing += sample_weight * treatments.row(sort_index).square();
      num_small_w_missing += (treatments.row(sort_index).transpose() < mean_node_w).cast<int>();
    } else {
      weight_sums[split_index] += sample_weight;
      sums.row(split_index) += sample_weight * responses_by_sample.row(sample);
      ++counter[split_index];

      sums_w.row(split_index) += sample_weight * treatments.row(sort_index);
      sums_w_squared.row(split_index) += sample_weight * treatments.row(sort_index).square();
      num_small_w.row(split_index) += (treatments.row(sort_index).transpose() < mean_node_w).cast<int>();
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
  Eigen::Ref<Eigen::ArrayXd> sum_left = sum_missing;
  Eigen::Ref<Eigen::ArrayXd> sum_left_w = sum_w_missing;
  Eigen::Ref<Eigen::ArrayXd> sum_left_w_squared = sum_w_squared_missing;
  Eigen::Ref<Eigen::ArrayXi> num_left_small_w = num_small_w_missing;

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
      sum_left.setZero();
      sum_left_w.setZero();
      sum_left_w_squared.setZero();
      num_left_small_w.setZero();
    }

    for (size_t i = 0; i < num_splits; ++i) {
      // not necessary to evaluate sending right when splitting on NaN.
      if (i == 0 && !send_left) {
        continue;
      }

      n_left += counter[i];
      weight_sum_left += weight_sums[i];
      num_left_small_w += num_small_w.row(i);
      sum_left += sums.row(i);
      sum_left_w += sums_w.row(i);
      sum_left_w_squared += sums_w_squared.row(i);

      // Skip this split if the left child does not contain enough
      // w values below and above the parent's mean.
      // We have Eigen::ArrayXi num_left_large_w = n_left - num_left_small_w; but write down the expressions
      // in-place so Eigen can resolve them with the minimum amount of unnecessary copies at compile time.
      // Same as: if ((num_left_small_w < min_node_size).any() || (num_left_large_w < min_node_size).any()
      if ((num_left_small_w < min_node_size).any() || (n_left - num_left_small_w < min_node_size).any()) {
        continue;
      }

      // Stop if the right child does not contain enough w values below
      // and above the parent's mean.
      size_t n_right = num_samples - n_left;
      // We have:
      // Eigen::ArrayXi num_right_small_w = num_node_small_w - num_left_small_w
      // Eigen::ArrayXi num_right_large_w = n_right - num_right_small_w
      // Same as: if ((num_right_small_w < min_node_size).any() || (num_right_large_w < min_node_size).any())
      if ((num_node_small_w - num_left_small_w < min_node_size).any() ||
          (n_right - num_node_small_w + num_left_small_w < min_node_size).any()) {
        break;
      }

      // Calculate relevant quantities for the left child.
      // We have:
      // Eigen::ArrayXd size_left = sum_left_w_squared - sum_left_w.square() / weight_sum_left
      // Skip this split if the left child's variance is too small.
      // Same as: if ((size_left < min_child_size).any() || (imbalance_penalty > 0.0 && (size_left == 0).all()))
      if ((sum_left_w_squared - sum_left_w.square() / weight_sum_left < min_child_size).any() ||
          (imbalance_penalty > 0.0 && (sum_left_w_squared - sum_left_w.square() / weight_sum_left == 0).all())) {
        continue;
      }

      // Calculate relevant quantities for the right child.
      double weight_sum_right = weight_sum_node - weight_sum_left;
      // We have:
      // Eigen::ArrayXd sum_right_w_squared = sum_node_w_squared - sum_left_w_squared;
      // Eigen::ArrayXd sum_right_w = sum_node_w - sum_left_w;
      // Eigen::ArrayXd size_right = sum_right_w_squared - sum_right_w.square() / weight_sum_right;
      // Skip this split if the right child's variance is too small.
      // Same as: if ((size_right < min_child_size).any() || (imbalance_penalty > 0.0 && (size_right == 0).all()))
      if ((sum_node_w_squared - sum_left_w_squared - (sum_node_w - sum_left_w).square() / weight_sum_right < min_child_size).any() ||
          (imbalance_penalty > 0.0 && (sum_node_w_squared - sum_left_w_squared - (sum_node_w - sum_left_w).square() / weight_sum_right == 0).all())) {
        continue;
      }

      // Calculate the decrease in impurity.
      // We have `Eigen::ArrayXd sum_right = sum_node - sum_left` but write down the expression
      // in-place below to avoid unnecessary temporaries.
      double decrease = sum_left.square().sum() / weight_sum_left +
                        (sum_node - sum_left).square().sum() / weight_sum_right;

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

} // namespace grf
