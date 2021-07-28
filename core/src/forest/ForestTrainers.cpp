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
#include "prediction/CausalSurvivalPredictionStrategy.h"
#include "prediction/InstrumentalPredictionStrategy.h"
#include "prediction/MultiCausalPredictionStrategy.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "prediction/MultiRegressionPredictionStrategy.h"
#include "prediction/ProbabilityPredictionStrategy.h"
#include "relabeling/CausalSurvivalRelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "relabeling/MultiCausalRelabelingStrategy.h"
#include "relabeling/LLRegressionRelabelingStrategy.h"
#include "relabeling/NoopRelabelingStrategy.h"
#include "relabeling/MultiNoopRelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "splitting/factory/InstrumentalSplittingRuleFactory.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "splitting/factory/RegressionSplittingRuleFactory.h"
#include "splitting/factory/MultiCausalSplittingRuleFactory.h"
#include "splitting/factory/MultiRegressionSplittingRuleFactory.h"
#include "splitting/factory/SurvivalSplittingRuleFactory.h"
#include "splitting/factory/CausalSurvivalSplittingRuleFactory.h"

namespace grf {

ForestTrainer instrumental_trainer(double reduced_form_weight,
                                   bool stabilize_splits) {

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(reduced_form_weight));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory = stabilize_splits
          ? std::unique_ptr<SplittingRuleFactory>(new InstrumentalSplittingRuleFactory())
          : std::unique_ptr<SplittingRuleFactory>(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());

  std::string verbose_operation_name = "Growing instrumental forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
}

ForestTrainer multi_causal_trainer(size_t num_treatments,
                                   size_t num_outcomes,
                                   bool stabilize_splits) {
  size_t response_length = num_treatments * num_outcomes;
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new MultiCausalRelabelingStrategy(response_length));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory = stabilize_splits
    ? std::unique_ptr<SplittingRuleFactory>(new MultiCausalSplittingRuleFactory(response_length, num_treatments))
    : std::unique_ptr<SplittingRuleFactory>(new MultiRegressionSplittingRuleFactory(response_length));
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new MultiCausalPredictionStrategy(num_treatments, num_outcomes));

  std::string verbose_operation_name = "Growing multi-causal forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
}

ForestTrainer quantile_trainer(const std::vector<double>& quantiles) {
    std::unique_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(
      new ProbabilitySplittingRuleFactory(quantiles.size() + 1));

  std::string verbose_operation_name = "Growing quantile forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       nullptr,
                       verbose_operation_name);
}

ForestTrainer probability_trainer(size_t num_classes) {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new ProbabilitySplittingRuleFactory(num_classes));
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new ProbabilityPredictionStrategy(num_classes));

  std::string verbose_operation_name = "Growing probability forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
}

ForestTrainer regression_trainer() {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  std::string verbose_operation_name = "Growing regression forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
}

ForestTrainer multi_regression_trainer(size_t num_outcomes) {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new MultiNoopRelabelingStrategy(num_outcomes));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new MultiRegressionSplittingRuleFactory(num_outcomes));
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new MultiRegressionPredictionStrategy(num_outcomes));

  std::string verbose_operation_name = "Growing multi-regression forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
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

  std::string verbose_operation_name = "Growing local linear regression forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
}

ForestTrainer survival_trainer() {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new SurvivalSplittingRuleFactory());

  std::string verbose_operation_name = "Growing survival forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       nullptr,
                       verbose_operation_name);
}

ForestTrainer causal_survival_trainer(bool stabilize_splits) {

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new CausalSurvivalRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory = stabilize_splits
          ? std::unique_ptr<SplittingRuleFactory>(new CausalSurvivalSplittingRuleFactory())
          : std::unique_ptr<SplittingRuleFactory>(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new CausalSurvivalPredictionStrategy());
  std::string verbose_operation_name = "Growing causal survival forest";

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy),
                       verbose_operation_name);
}

} // namespace grf
