/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

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

#include <unordered_map>
#include "SplittingRule.h"
#include "ProbabilitySplittingRuleFactory.h"
#include "ProbabilitySplittingRule.h"


ProbabilitySplittingRuleFactory::ProbabilitySplittingRuleFactory(Data *data, size_t num_classes):
    data(data), num_classes(num_classes) {}

std::shared_ptr<SplittingRule> ProbabilitySplittingRuleFactory::create() {
  return std::shared_ptr<SplittingRule>(new ProbabilitySplittingRule(data, num_classes));
}
