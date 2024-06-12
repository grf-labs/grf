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

#include "splitting/factory/MultiCausalSplittingRuleFactory.h"
#include "splitting/MultiCausalSplittingRule.h"

namespace grf {

MultiCausalSplittingRuleFactory::MultiCausalSplittingRuleFactory(size_t response_length,
                                                                 size_t num_treatments):
  response_length(response_length),
  num_treatments(num_treatments) {}

std::unique_ptr<SplittingRule> MultiCausalSplittingRuleFactory::create(size_t max_num_unique_values,
                                                                       const TreeOptions& options) const {
  return std::unique_ptr<SplittingRule>(new MultiCausalSplittingRule(
      max_num_unique_values,
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      response_length,
      num_treatments));
}

} // namespace grf
