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

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
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

ForestTrainer probability_trainer(size_t num_classes) {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new ProbabilitySplittingRuleFactory(num_classes));
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new ProbabilityPredictionStrategy(num_classes));

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
}

ForestTrainer regression_trainer() {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
}

ForestTrainer multi_regression_trainer(size_t num_outcomes) {
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new MultiNoopRelabelingStrategy(num_outcomes));
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory(new MultiRegressionSplittingRuleFactory(num_outcomes));
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new MultiRegressionPredictionStrategy(num_outcomes));

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

ForestTrainer causal_survival_trainer(bool stabilize_splits) {

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new CausalSurvivalRelabelingStrategy());
  std::unique_ptr<SplittingRuleFactory> splitting_rule_factory = stabilize_splits
          ? std::unique_ptr<SplittingRuleFactory>(new CausalSurvivalSplittingRuleFactory())
          : std::unique_ptr<SplittingRuleFactory>(new RegressionSplittingRuleFactory());
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new CausalSurvivalPredictionStrategy());

  return ForestTrainer(std::move(relabeling_strategy),
                       std::move(splitting_rule_factory),
                       std::move(prediction_strategy));
}

} // namespace grf
