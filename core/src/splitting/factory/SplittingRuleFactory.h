#ifndef GRADIENTFOREST_SPLITTINGRULEFACTORY_H
#define GRADIENTFOREST_SPLITTINGRULEFACTORY_H

#include "Data.h"
#include "SplittingRule.h"

class SplittingRuleFactory {
public:
  virtual SplittingRule *create() = 0;
};

#endif //GRADIENTFOREST_SPLITTINGRULEFACTORY_H
