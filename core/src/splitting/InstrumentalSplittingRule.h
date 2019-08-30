/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#ifndef GRF_INSTRUMENTALSPLITTINGRULE_H
#define GRF_INSTRUMENTALSPLITTINGRULE_H

#include "commons/Data.h"
#include "commons/Data.h"
#include "splitting/SplittingRule.h"

class InstrumentalSplittingRule final: public SplittingRule {
public:
  InstrumentalSplittingRule(const Data* data,
                            uint min_node_size,
                            double alpha,
                            double imbalance_penalty);
  ~InstrumentalSplittingRule();

  bool find_best_split(size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const std::unordered_map<size_t, double>& labels_by_sample,
                       const std::vector<std::vector<size_t>>& samples,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values);

private:
  void find_best_split_value_small_q(size_t node,
                                     size_t var,
                                     size_t num_samples,
                                     double sum_node,
                                     double mean_node_z,
                                     size_t num_node_small_z,
                                     double sum_node_w,
                                     double sum_node_w_squared,
                                     double min_child_size,
                                     double& best_value,
                                     size_t& best_var,
                                     double& best_decrease,
                                     const std::unordered_map<size_t, double>& responses_by_sample,
                                     const std::vector<std::vector<size_t>>& samples);
  void find_best_split_value_large_q(size_t node,
                                     size_t var,
                                     size_t num_samples,
                                     double sum_node,
                                     double mean_node_z,
                                     size_t num_node_small_z,
                                     double sum_node_w,
                                     double sum_node_w_squared,
                                     double min_child_size,
                                     double& best_value,
                                     size_t& best_var,
                                     double& best_decrease,
                                     const std::unordered_map<size_t, double>& responses_by_sample,
                                     const std::vector<std::vector<size_t>>& samples);

  const Data* data;

  size_t* counter;
  double* sums;
  size_t* num_small_z;
  double* sums_z;
  double* sums_z_squared;

  uint min_node_size;
  double alpha;
  double imbalance_penalty;

  DISALLOW_COPY_AND_ASSIGN(InstrumentalSplittingRule);
};


#endif //GRF_INSTRUMENTALSPLITTINGRULE_H
