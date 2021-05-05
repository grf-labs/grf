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

#ifndef GRF_MULTICAUSALSPLITTINGRULE_H
#define GRF_MULTICAUSALSPLITTINGRULE_H

#include "commons/Data.h"
#include "splitting/SplittingRule.h"

namespace grf {

/**
 * This splitting rule is a multivariate extension of {@link InstrumentalSplittingRule}.
 * The contraints are
 * (1)
 * Define A_ik = W_ik < \bar{W}_k,
 * where k = 1,...,K (number of treatments) and \bar(W)_k is the average of the k'th
 * treatment in the parent.
 * Requirement: for each k, each child node must contain at least min_node_size
 * observations with A_ik = 1, and A_ik = 0.
 *
 * (2)
 * Define size(parent_k) = \sum_i{i in node} (W_ik - \bar{W}_k)^2.
 * Requirement: for each k, the size of each child must be greater or equal to
 * size(parent_k) * alpha, where alpha: (0, 0.25) is a tuning parameter.
 *
 */
class MultiCausalSplittingRule final: public SplittingRule {
public:
  MultiCausalSplittingRule(size_t max_num_unique_values,
                           uint min_node_size,
                           double alpha,
                           double imbalance_penalty,
                           size_t response_length,
                           size_t num_treatments);

  ~MultiCausalSplittingRule();

  bool find_best_split(const Data& data,
                       size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const Eigen::ArrayXXd& responses_by_sample,
                       const std::vector<std::vector<size_t>>& samples,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values,
                       std::vector<bool>& send_missing_left);

private:
  void find_best_split_value(const Data& data,
                             size_t node,
                             size_t var,
                             size_t num_samples,
                             double weight_sum_node,
                             const Eigen::ArrayXd& sum_node,
                             const Eigen::ArrayXd& mean_node_w,
                             const Eigen::ArrayXi& sum_node_small_w,
                             const Eigen::ArrayXd& sum_node_w,
                             const Eigen::ArrayXd& sum_node_w_squared,
                             const Eigen::ArrayXd& min_child_size,
                             const Eigen::ArrayXXd& treatments,
                             double& best_value,
                             size_t& best_var,
                             double& best_decrease,
                             bool& best_send_missing_left,
                             const Eigen::ArrayXXd& responses_by_sample,
                             const std::vector<std::vector<size_t>>& samples);

  size_t* counter;
  double* weight_sums;
  Eigen::ArrayXXd sums;
  Eigen::ArrayXXi num_small_w;
  Eigen::ArrayXXd sums_w;
  Eigen::ArrayXXd sums_w_squared;

  uint min_node_size;
  double alpha;
  double imbalance_penalty;
  size_t response_length;
  size_t num_treatments;

  DISALLOW_COPY_AND_ASSIGN(MultiCausalSplittingRule);
};

} // namespace grf

#endif //GRF_MULTICAUSALSPLITTINGRULE_H
