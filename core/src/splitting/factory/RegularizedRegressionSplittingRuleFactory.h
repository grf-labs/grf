/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#ifndef GRF_REGULARIZEDREGRESSIONSPLITTINGRULEFACTORY_H
#define GRF_REGULARIZEDREGRESSIONSPLITTINGRULEFACTORY_H


#include "commons/DefaultData.h"
#include "splitting/SplittingRule.h"
#include "SplittingRuleFactory.h"

class RegularizedRegressionSplittingRuleFactory: public SplittingRuleFactory {
public:
  /**
   * Creates a factory that produces regularized regression splitting rules.
   *
   * This splitting rule applies a penalty to avoid splits too close to the
   * edge of the node's data.
   *
   * data: A pointer to the training data.
   * lambda: A tuning parameter to control the amount of regularization applied
   *     to each split.
   * downweight_penalty: Whether or not to downweight the penalty (experimental,
   *     will be removed once we run tuning experiments).
   */
  RegularizedRegressionSplittingRuleFactory(Data* data, double lambda, bool downweight_penalty);
  std::shared_ptr<SplittingRule> create();

private:
  Data* data;
  double lambda;
  bool downweight_penalty;
};


#endif //GRF_REGULARIZEDREGRESSIONSPLITTINGRULEFACTORY_H
