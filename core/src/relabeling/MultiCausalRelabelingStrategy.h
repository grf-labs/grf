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

#ifndef GRF_MULTICAUSALRELABELINGSTRATEGY_H
#define GRF_MULTICAUSALRELABELINGSTRATEGY_H

#include <vector>

#include "commons/Data.h"
#include "relabeling/RelabelingStrategy.h"
#include "tree/Tree.h"

namespace grf {

/**
 * This relabeling strategy is a multi-treatment extension of {@link InstrumentalRelabelingStrategy}.
 * We compute the vector-valued gradient for tau = [tau_1, ..., tau_K] wrt. observation i in the regression
 * Y = c + tau W + e,
 * where W is a vector of treatment variables (k = 1,...,K). See equation (20) in
 * https://arxiv.org/pdf/1610.01271.pdf. The response Y may be multivariate (m = 1,..,M),
 * in which case we concatenate the influence vectors for each response.
 *
 * The output of this method is a vector with coefficient influences for each observation ordered according to:
 * [\delta \tau_{11}, ..., \delta \tau_{1K}, ..., \delta \tau_{M1}, ..., \delta \tau_{MK}] * gradient_weight,
 * where `*` denotes elementwise multiplication with the optional weight vector `gradient_weight` (by default
 * equal to [1, ..., 1]).
 *
 */
class MultiCausalRelabelingStrategy final: public RelabelingStrategy {
public:
  MultiCausalRelabelingStrategy(size_t response_length,
                                const std::vector<double>& gradient_weights);

  bool relabel(
      const std::vector<size_t>& samples,
      const Data& data,
      Eigen::ArrayXXd& responses_by_sample) const;

  size_t get_response_length() const;

private:
  size_t response_length;
  std::vector<double> gradient_weights;
};

} // namespace grf

#endif //GRF_MULTICAUSALRELABELINGSTRATEGY_H
