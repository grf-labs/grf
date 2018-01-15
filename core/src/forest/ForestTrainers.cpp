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
#include "prediction/CustomPredictionStrategy.h"
#include "prediction/InstrumentalPredictionStrategy.h"
#include "prediction/QuantilePredictionStrategy.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "relabeling/CustomRelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "relabeling/NoopRelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "splitting/factory/RegressionSplittingRuleFactory.h"
#include "splitting/factory/RegularizedRegressionSplittingRuleFactory.h"

ForestTrainer ForestTrainers::instrumental_trainer(size_t outcome_index,
                                                   size_t treatment_index,
                                                   size_t instrument_index,
                                                   double split_regularization,
                                                   double alpha,
                                                   const ForestOptions& options) {
  std::unordered_map<size_t, size_t> observables = {
      {Observations::OUTCOME, outcome_index},
      {Observations::TREATMENT, treatment_index},
      {Observations::INSTRUMENT, instrument_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(split_regularization));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(alpha));
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());

  return ForestTrainer(observables, options, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::quantile_trainer(size_t outcome_index,
                                               const std::vector<double>& quantiles,
                                               double alpha,
                                               const ForestOptions& options) {
  std::unordered_map<size_t, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(
      new ProbabilitySplittingRuleFactory(alpha, quantiles.size() + 1));

  return ForestTrainer(observables, options, relabeling_strategy, splitting_rule_factory, NULL);
}

ForestTrainer ForestTrainers::regression_trainer(size_t outcome_index,
                                                 double alpha,
                                                 const ForestOptions& options) {
  std::unordered_map<size_t, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(alpha));
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(observables, options, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::custom_trainer(size_t outcome_index,
                                             double alpha,
                                             const ForestOptions& options) {
  std::unordered_map<size_t, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new CustomRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(alpha));

  return ForestTrainer(observables, options, relabeling_strategy, splitting_rule_factory, NULL);
}

ForestTrainer ForestTrainers::regularized_regression_trainer(size_t outcome_index,
                                                             double lambda,
                                                             bool downweight_penalty,
                                                             const ForestOptions& options) {
  std::unordered_map<size_t, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(
      new RegularizedRegressionSplittingRuleFactory(lambda, downweight_penalty));
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(observables, options, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::regularized_instrumental_trainer(size_t outcome_index,
                                                               size_t treatment_index,
                                                               size_t instrument_index,
                                                               double split_regularization,
                                                               double lambda,
                                                               bool downweight_penalty,
                                                               const ForestOptions& options) {
    std::unordered_map<size_t, size_t> observables = {
      {Observations::OUTCOME, outcome_index},
      {Observations::TREATMENT, treatment_index},
      {Observations::INSTRUMENT, instrument_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(split_regularization));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(
      new RegularizedRegressionSplittingRuleFactory(lambda, downweight_penalty));
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());

  return ForestTrainer(observables, options, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

