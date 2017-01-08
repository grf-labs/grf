#include <unordered_map>
#include "SplittingRule.h"
#include "ProbabilitySplittingRuleFactory.h"
#include "ProbabilitySplittingRule.h"


ProbabilitySplittingRuleFactory::ProbabilitySplittingRuleFactory(Data *data, size_t num_classes):
    data(data), num_classes(num_classes) {}

std::shared_ptr<SplittingRule> ProbabilitySplittingRuleFactory::create() {
  return std::shared_ptr<SplittingRule>(new ProbabilitySplittingRule(data, num_classes));
}
