#include "ForestTrainers.h"
#include "InstrumentalRelabelingStrategy.h"
#include "RegressionSplittingRuleFactory.h"
#include "ProbabilitySplittingRuleFactory.h"
#include "QuantileRelabelingStrategy.h"
#include "NoopRelabelingStrategy.h"

ForestTrainer ForestTrainers::instrumental_trainer(Data* data,
                                                   size_t outcome_index,
                                                   size_t treatment_index,
                                                   size_t instrument_index) {
  std::unordered_map<std::string, size_t> observables = {
      {Observations::OUTCOME, outcome_index},
      {Observations::TREATMENT, treatment_index},
      {Observations::INSTRUMENT, instrument_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(data));

  return ForestTrainer(observables, relabeling_strategy, splitting_rule_factory);
}

ForestTrainer ForestTrainers::quantile_trainer(Data* data,
                                               size_t outcome_index,
                                               const std::vector<double>& quantiles) {
  std::unordered_map<std::string, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(
      new ProbabilitySplittingRuleFactory(data, quantiles.size() + 1));

  return ForestTrainer(observables, relabeling_strategy, splitting_rule_factory);
}

ForestTrainer ForestTrainers::regression_trainer(Data* data,
                                                 size_t outcome_index) {
  std::unordered_map<std::string, size_t> observables = {{Observations::OUTCOME, outcome_index}};

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory(new RegressionSplittingRuleFactory(data));

  return ForestTrainer(observables, relabeling_strategy, splitting_rule_factory);
}
