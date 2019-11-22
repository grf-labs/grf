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

#include "commons/DefaultData.h"
#include "splitting/SplittingRule.h"
#include "tree/Tree.h"

namespace grf {

class RegressionSplittingRule final: public SplittingRule {
public:
  RegressionSplittingRule(size_t max_num_unique_values,
                          double alpha,
                          double imbalance_penalty);

  ~RegressionSplittingRule();

  /**
   * Finds the best split at a given node in the tree.
   *
   * Is called repeatedly to build a tree in a breadth-first fashion.
   *
   * @param data: the data matrix containing all test samples.
   * @param node: the node id in the tree.
   * @param possible_split_vars: a vector of valid covariate IDs.
   * @param responses_by_sample: a map from sample ID to response.
   * @param samples: a vector of samples at the given node.
   * @param split_vars: the output of the method, the best split variable, stored at node.
   * @param split_values: the output of the method, the best split value, stored at node.
   * @return a boolean that will be true if no best split was found.
   *
   * Details:
   *
   * At each split variable j, this method calls into `find_best_split_value_small_q` or
   * `find_best_split_value_large_q` depending on the following ratio q:
   * the number of samples in the node (nj) over the total number of unique values in the data at variable j (Nj).
   * If this value is less than Q_THRESHOLD (0.02 by default) then the small_q method is called,
   * otherwise large_q.
   *
   * An expensive computation in finding the best split is sorting all the values in order to place samples
   * on the left or right of the split.
   * `small_q` sorts the data in the node when it is called, which has time complexity O(nj log nj), then
   * iterates over all nj points calculating the decrease in impurity.
   * `large_q` accesses a global sort order for all Nj points created at the initalization of forest training.
   * To find the best split at the nj points in the node, a full scan
   * over all Nj points is done, getting the proper position of the sample. This is dominated by O(Nj).
   *
   * The two splitting strategies balances O(nj log nj) and O(Nj). In large nodes at the root of the tree,
   * Nj is smaller than nj log nj, but in smaller nodes, deeper down the tree, nj is small and
   * Nj is larger than nj log nj.
   *
   * The exact value of Q_THRESHOLD has been determined empirically by the ranger developers.
   */
  bool find_best_split(const Data& data,
                       size_t node,
                       const std::vector<size_t>& possible_split_vars,
                       const std::vector<double>& responses_by_sample,
                       const std::vector<std::vector<size_t>>& samples,
                       std::vector<size_t>& split_vars,
                       std::vector<double>& split_values);

private:
  void find_best_split_value_small_q(const Data& data,
                                     size_t node,
                                     size_t var,
                                     double sum_node,
                                     size_t size_node,
                                     size_t min_child_size,
                                     double& best_value,
                                     size_t& best_var,
                                     double& best_decrease,
                                     const std::vector<double>& responses_by_sample,
                                     const std::vector<std::vector<size_t>>& samples);
  void find_best_split_value_large_q(const Data& data,
                                     size_t node,
                                     size_t var,
                                     double sum_node,
                                     size_t size_node,
                                     size_t mind_child_size,
                                     double& best_value,
                                     size_t& best_var,
                                     double& best_decrease,
                                     const std::vector<double>& responses_by_sample,
                                     const std::vector<std::vector<size_t>>& samples);

  size_t* counter;
  double* sums;

  double alpha;
  double imbalance_penalty;

  DISALLOW_COPY_AND_ASSIGN(RegressionSplittingRule);
};

} // namespace grf

#endif //GRF_REGRESSIONSPLITTINGRULE_H
