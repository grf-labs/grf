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

#include "forest/ForestTrainers.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "prediction/InstrumentalPredictionStrategy.h"
#include "splitting/factory/RegressionSplittingRuleFactory.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "prediction/QuantilePredictionStrategy.h"
#include "relabeling/NoopRelabelingStrategy.h"

ForestTrainer ForestTrainers::instrumental_trainer(Data* data,
                                                   size_t outcome_index,
                                                   size_t treatment_index,
                                                   size_t instrument_index,
                                                   double split_regularization) {
  std::unordered_map<size_t, size_t> observables = {
      {Observations::OUTCOME, outcome_index},
      {Observations::TREATMENT, treatment_index},
      {Observations::INSTRUMENT, instrument_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(split_regularization));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(data));
  std::shared_ptr<PredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());

  return ForestTrainer(observables, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::quantile_trainer(Data* data,
                                               size_t outcome_index,
                                               const std::vector<double>& quantiles) {
  std::unordered_map<size_t, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(
      new ProbabilitySplittingRuleFactory(data, quantiles.size() + 1));
  std::shared_ptr<PredictionStrategy> prediction_strategy(new QuantilePredictionStrategy(quantiles));

  return ForestTrainer(observables, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}

ForestTrainer ForestTrainers::regression_trainer(Data* data,
                                                 size_t outcome_index) {
  std::unordered_map<size_t, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(data));
  std::shared_ptr<PredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());

  return ForestTrainer(observables, relabeling_strategy, splitting_rule_factory, prediction_strategy);
}
