#ifndef GRADIENTFOREST_SPLITTINGRULEFACTORY_H
#define GRADIENTFOREST_SPLITTINGRULEFACTORY_H

#include <memory>

#include "Data.h"
#include "SplittingRule.h"

class SplittingRuleFactory {
public:
  virtual std::shared_ptr<SplittingRule> create() = 0;
};

#endif //GRADIENTFOREST_SPLITTINGRULEFACTORY_H
