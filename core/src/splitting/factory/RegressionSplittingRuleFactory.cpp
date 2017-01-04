#include "RegressionSplittingRuleFactory.h"
#include "RegressionSplittingRule.h"

RegressionSplittingRuleFactory::RegressionSplittingRuleFactory(Data *data):
    data(data) {}

SplittingRule *RegressionSplittingRuleFactory::create() {
  return new RegressionSplittingRule(data);
}
