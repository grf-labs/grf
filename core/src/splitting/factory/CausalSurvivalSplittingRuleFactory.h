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

#ifndef GRF_CAUSALSURVIVALSPLITTINGRULEFACTORY_H
#define GRF_CAUSALSURVIVALSPLITTINGRULEFACTORY_H


#include "splitting/factory/SplittingRuleFactory.h"

namespace grf {

class CausalSurvivalSplittingRuleFactory final: public SplittingRuleFactory {
public:
  CausalSurvivalSplittingRuleFactory() = default;
  std::unique_ptr<SplittingRule> create(size_t max_num_unique_values,
                                        const TreeOptions& options) const;
private:
  DISALLOW_COPY_AND_ASSIGN(CausalSurvivalSplittingRuleFactory);
};

} // namespace grf

#endif //GRF_CAUSALSURVIVALSPLITTINGRULEFACTORY_H
