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

#ifndef GRF_PROBABILITYSPLITTINGRULEFACTORY_H
#define GRF_PROBABILITYSPLITTINGRULEFACTORY_H

#include <vector>

#include "commons/globals.h"
#include "commons/DefaultData.h"
#include "splitting/factory/SplittingRuleFactory.h"

class ProbabilitySplittingRuleFactory: public SplittingRuleFactory {
public:
  /**
   * Creates a factory that produces standard classification splitting rules.
   *
   * data: A pointer to the training data.
   * alpha: The minimum fraction of samples that are allowed to be on either
   *     side of the split. Splits that are too uneven according to this
   *     parameter will not be considered.
   * num_classes: The number of unique classes in the data.
   */
  ProbabilitySplittingRuleFactory(double alpha, double lambda, size_t num_classes);
  std::shared_ptr<SplittingRule> create(Data* data);

private:
  double alpha;
  double lambda;
  size_t num_classes;

  DISALLOW_COPY_AND_ASSIGN(ProbabilitySplittingRuleFactory);
};

#endif //GRF_PROBABILITYSPLITTINGRULEFACTORY_H
