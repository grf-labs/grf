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

#ifndef GRF_SURVIVALSPLITTINGRULE_H
#define GRF_SURVIVALSPLITTINGRULE_H

#include "Eigen/Dense"

#include "commons/Data.h"
#include "splitting/SplittingRule.h"
#include "tree/Tree.h"

namespace grf {

class SurvivalSplittingRule final: public SplittingRule {
public:
  SurvivalSplittingRule(double alpha);

  bool find_best_split(const Data& data,
                       size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const Eigen::ArrayXXd& responses_by_sample,
                       const std::vector<std::vector<size_t>>& samples_by_node,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values,
                       std::vector<bool>& send_missing_left);

 /**
  * This member is public for unit testing purposes. It returns an additional
  * output value, the best logrank statistic.
  */
 void find_best_split_internal(const Data& data,
                               const std::vector<size_t>& possible_split_vars,
                               const Eigen::ArrayXXd& responses_by_sample,
                               const std::vector<size_t>& samples,
                               double& best_value,
                               size_t& best_var,
                               bool& best_send_missing_left,
                               double& best_logrank);

private:
  void find_best_split_value(const Data& data,
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
                             const std::vector<double>& count_failure,
                             const std::vector<double>& at_risk,
                             const std::vector<double>& numerator_weights,
                             const std::vector<double>& denominator_weights);

  inline double compute_logrank(size_t num_failures,
                                size_t n_left,
                                std::vector<double>& cum_sums,
                                const std::vector<double>& left_count_failure,
                                const std::vector<double>& left_count_censor,
                                const std::vector<double>& at_risk,
                                const std::vector<double>& numerator_weights,
                                const std::vector<double>& denominator_weights);

  double alpha;

  DISALLOW_COPY_AND_ASSIGN(SurvivalSplittingRule);
};

} // namespace grf

#endif //GRF_SURVIVALSPLITTINGRULE_H
