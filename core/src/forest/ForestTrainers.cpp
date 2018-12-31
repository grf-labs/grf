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

#include "forest/ForestTrainers.h"
#include "prediction/InstrumentalPredictionStrategy.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "relabeling/CustomRelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "relabeling/NoopRelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "splitting/factory/InstrumentalSplittingRuleFactory.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "splitting/factory/RegressionSplittingRuleFactory.h"


ForestTrainer ForestTrainers::instrumental_trainer(double reduced_form_weight,
                                                   bool stabilize_splits) {

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(reduced_form_weight));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory = stabilize_splits
          ? std::shared_ptr<SplittingRuleFactory>(new InstrumentalSplittingRuleFactory())
          : std::shared_ptr<SplittingRuleFactory>(new RegressionSplittingRuleFactory());
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());

  return ForestTrainer(relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::quantile_trainer(const std::vector<double>& quantiles) {
  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(
      new ProbabilitySplittingRuleFactory(quantiles.size() + 1));

  return ForestTrainer(relabeling_strategy, splitting_rule_factory, NULL);
}

ForestTrainer ForestTrainers::regression_trainer() {
  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::custom_trainer() {
  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new CustomRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());

  return ForestTrainer(relabeling_strategy, splitting_rule_factory, NULL);
}
