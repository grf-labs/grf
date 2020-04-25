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

#ifndef GRF_SPLITTINGRULE_H
#define GRF_SPLITTINGRULE_H

#include <vector>

#include "commons/Data.h"

namespace grf {

class SplittingRule {
public:
  virtual ~SplittingRule() {}

  /**
   * Finds the best split at a given node in the tree.
   *
   * Is called repeatedly to build a tree in a breadth-first fashion.
   *
   * @param data: the data matrix containing all test samples.
   * @param node: the node id in the tree.
   * @param possible_split_vars: a vector of valid covariate IDs.
   * @param responses_by_sample: the response for each sample.
   * @param samples: a vector of samples at the given node.
   * @param split_vars: the output of the method, the best split variable, stored at node.
   * @param split_values: the output of the method, the best split value, stored at node.
   * @return a boolean that will be true if no best split was found.
   *
   */
  virtual bool find_best_split(const Data& data,
                               size_t node,
                               const std::vector<size_t>& possible_split_vars,
                               const std::vector<double>& responses_by_sample,
                               const std::vector<std::vector<size_t>>& samples,
                               std::vector<size_t>& split_vars,
                               std::vector<double>& split_values,
                               std::vector<bool>& send_missing_left) = 0;
};

} // namespace grf

#endif //GRF_SPLITTINGRULE_H
