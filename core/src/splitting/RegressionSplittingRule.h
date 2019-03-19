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

#ifndef GRF_REGRESSIONSPLITTINGRULE_H
#define GRF_REGRESSIONSPLITTINGRULE_H

#include "tree/Tree.h"
#include "splitting/SplittingRule.h"
#include <unordered_map>
#include "commons/DefaultData.h"

class RegressionSplittingRule: public SplittingRule {
public:
  RegressionSplittingRule(const Data* data,
                          double alpha,
                          double imbalance_penalty);

  ~RegressionSplittingRule();

  bool find_best_split(size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const std::unordered_map<size_t, double>& labels_by_sample,
                       const std::vector<std::vector<size_t>>& samples,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values);

private:
  void find_best_split_value_small_q(size_t node,
                                     size_t var,
                                     double sum_node,
                                     size_t size_node,
                                     size_t min_child_size,
                                     double& best_value,
                                     size_t& best_var,
                                     double& best_decrease,
                                     const std::unordered_map<size_t, double>& labels_by_sample,
                                     const std::vector<std::vector<size_t>>& samples);
  void find_best_split_value_large_q(size_t node,
                                     size_t var,
                                     double sum_node,
                                     size_t size_node,
                                     size_t mind_child_size,
                                     double& best_value,
                                     size_t& best_var,
                                     double& best_decrease,
                                     const std::unordered_map<size_t, double>& responses_by_sample,
                                     const std::vector<std::vector<size_t>>& samples);

  const Data* data;
  size_t* counter;
  double* sums;

  double alpha;
  double imbalance_penalty;

  DISALLOW_COPY_AND_ASSIGN(RegressionSplittingRule);
};


#endif //GRF_REGRESSIONSPLITTINGRULE_H
