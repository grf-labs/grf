#include <unordered_map>
#include <splitting/SplittingRule.h>
#include "ProbabilitySplittingRuleFactory.h"
#include "ProbabilitySplittingRule.h"


ProbabilitySplittingRuleFactory::ProbabilitySplittingRuleFactory(Data *data, size_t num_classes):
    data(data), num_classes(num_classes) {}

SplittingRule *ProbabilitySplittingRuleFactory::create() {
  return new ProbabilitySplittingRule(data, num_classes);
}
