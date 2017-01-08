#include "RegressionSplittingRuleFactory.h"
#include "RegressionSplittingRule.h"

RegressionSplittingRuleFactory::RegressionSplittingRuleFactory(Data *data):
    data(data) {}

std::shared_ptr<SplittingRule> RegressionSplittingRuleFactory::create() {
  return std::shared_ptr<SplittingRule>(new RegressionSplittingRule(data));
}
