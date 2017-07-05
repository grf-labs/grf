/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_REGULARIZEDREGRESSIONSPLITTINGRULE_H
#define GRF_REGULARIZEDREGRESSIONSPLITTINGRULE_H


#include "commons/Data.h"
#include "SplittingRule.h"
#include <unordered_map>
#include <vector>

class RegularizedRegressionSplittingRule: public SplittingRule {
public:
  RegularizedRegressionSplittingRule(Data* data, double lambda, bool downweight_penalty);

  ~RegularizedRegressionSplittingRule();

  bool find_best_split(size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const std::unordered_map<size_t, double>& labels_by_sample,
                       const std::vector<std::vector<size_t>>& samples,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values);

private:
  virtual void find_best_split_value_small_q(size_t node,
                                             size_t var,
                                             double sum_node,
                                             size_t num_samples_node,
                                             double node_impurity,
                                             double& best_value,
                                             size_t& best_var,
                                             double& best_decrease,
                                             const std::unordered_map<size_t, double>& responses_by_sample,
                                             const std::vector<std::vector<size_t>>& samples);
  virtual void find_best_split_value_large_q(size_t node,
                                             size_t var,
                                             double sum_node,
                                             size_t num_samples_node,
                                             double node_impurity,
                                             double& best_value,
                                             size_t& best_var,
                                             double& best_decrease,
                                             const std::unordered_map<size_t, double>& responses_by_sample,
                                             const std::vector<std::vector<size_t>>& samples);

  Data* data;
  size_t* counter;
  double* sums;
  double lambda;
  bool downweight_penalty;

  DISALLOW_COPY_AND_ASSIGN(RegularizedRegressionSplittingRule);
};


#endif //GRF_REGULARIZEDREGRESSIONSPLITTINGRULE_H
