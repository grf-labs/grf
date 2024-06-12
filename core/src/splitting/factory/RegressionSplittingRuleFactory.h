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

#ifndef GRF_REGRESSIONSPLITTINGRULEFACTORY_H
#define GRF_REGRESSIONSPLITTINGRULEFACTORY_H


#include "splitting/factory/SplittingRuleFactory.h"

namespace grf {

/**
 * A factory that produces standard regression splitting rules.
 *
 * In addition to performing standard regression splits, this rule applies
 * a penalty to avoid splits too close to the edge of the node's data.
 */
class RegressionSplittingRuleFactory final: public SplittingRuleFactory {
public:
  RegressionSplittingRuleFactory() = default;
  std::unique_ptr<SplittingRule> create(size_t max_num_unique_values,
                                        const TreeOptions& options) const;
private:
  DISALLOW_COPY_AND_ASSIGN(RegressionSplittingRuleFactory);
};

} // namespace grf

#endif //GRF_REGRESSIONSPLITTINGRULEFACTORY_H
