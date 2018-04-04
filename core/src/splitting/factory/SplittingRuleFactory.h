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

#ifndef GRF_SPLITTINGRULEFACTORY_H
#define GRF_SPLITTINGRULEFACTORY_H

#include <memory>

#include "commons/Data.h"
#include "commons/Observations.h"
#include "splitting/SplittingRule.h"
#include "tree/TreeOptions.h"


class SplittingRuleFactory {
public:
  virtual std::shared_ptr<SplittingRule> create(Data* data,
                                                const Observations& observations,
                                                const TreeOptions& options) = 0;
};

#endif //GRF_SPLITTINGRULEFACTORY_H
