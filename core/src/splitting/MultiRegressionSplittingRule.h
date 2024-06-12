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

#ifndef GRF_MULTIREGRESSIONSPLITTINGRULE_H
#define GRF_MULTIREGRESSIONSPLITTINGRULE_H

#include "commons/Data.h"
#include "splitting/SplittingRule.h"
#include "tree/Tree.h"

namespace grf {

/**
 * This splitting rule is identical to `RegressionSplittingRule` with the extension
 * that it allows for splits on vector valued Y (responses_by_sample) by considering
 * the squared L2 norm of the response sums.
 */
class MultiRegressionSplittingRule final: public SplittingRule {
public:
  MultiRegressionSplittingRule(size_t max_num_unique_values,
                               double alpha,
                               double imbalance_penalty,
                               size_t num_outcomes);

  ~MultiRegressionSplittingRule();

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
                             double weight_sum_node,
                             const Eigen::ArrayXd& sum_node,
                             size_t size_node,
                             size_t min_child_size,
                             double& best_value,
                             size_t& best_var,
                             double& best_decrease,
                             bool& best_send_missing_left,
                             const Eigen::ArrayXXd& responses_by_sample,
                             const std::vector<std::vector<size_t>>& samples);

  size_t* counter;
  Eigen::ArrayXXd sums;
  double* weight_sums;

  double alpha;
  double imbalance_penalty;
  size_t num_outcomes;

  DISALLOW_COPY_AND_ASSIGN(MultiRegressionSplittingRule);
};

} // namespace grf

#endif //GRF_MULTIREGRESSIONSPLITTINGRULE_H
