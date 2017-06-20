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

#include <unordered_map>
#include "splitting/SplittingRule.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "splitting/ProbabilitySplittingRule.h"


ProbabilitySplittingRuleFactory::ProbabilitySplittingRuleFactory(Data* data,
                                                                 double alpha,
                                                                 size_t num_classes):
    data(data), alpha(alpha), num_classes(num_classes) {}

std::shared_ptr<SplittingRule> ProbabilitySplittingRuleFactory::create() {
  return std::shared_ptr<SplittingRule>(new ProbabilitySplittingRule(data, alpha, num_classes));
}
