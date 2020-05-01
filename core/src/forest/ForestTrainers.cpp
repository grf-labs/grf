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
#include "relabeling/LLRegressionRelabelingStrategy.h"
#include "relabeling/NoopRelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "splitting/factory/InstrumentalSplittingRuleFactory.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "splitting/factory/RegressionSplittingRuleFactory.h"
#include "splitting/factory/SurvivalSplittingRuleFactory.h"

namespace grf {

ForestTrainer instrumental_trainer(double reduced_form_weight,
                                   bool stabilize_splits) {

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(reduced_form_weight));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory = stabilize_splits
          ? std::unique_ptr<SplittingRuleFactory>(new InstrumentalSplittingRuleFactory())
          : std::unique_ptr<SplittingRuleFactory>(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
}

ForestTrainer quantile_trainer(const std::vector<double>& quantiles) {
    std::unique_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(
      new ProbabilitySplittingRuleFactory(quantiles.size() + 1));

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       nullptr);
}

ForestTrainer regression_trainer() {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
}

ForestTrainer ll_regression_trainer(double split_lambda,
                                   bool weight_penalty,
                                   const std::vector<double>& overall_beta,
                                   size_t ll_split_cutoff,
                                   std::vector<size_t> ll_split_variables) {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new LLRegressionRelabelingStrategy(split_lambda, weight_penalty, overall_beta,
                                                                                             ll_split_cutoff, ll_split_variables));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
}

ForestTrainer survival_trainer() {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new SurvivalSplittingRuleFactory());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       nullptr);
}

ForestTrainer custom_trainer() {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new CustomRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       nullptr);
}

} // namespace grf
