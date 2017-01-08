#ifndef GRADIENTFOREST_REGRESSIONSPLITTINGRULEFACTORY_H
#define GRADIENTFOREST_REGRESSIONSPLITTINGRULEFACTORY_H


#include "SplittingRuleFactory.h"

class RegressionSplittingRuleFactory: public SplittingRuleFactory {
public:
  RegressionSplittingRuleFactory(Data* data);
  std::shared_ptr<SplittingRule> create();

private:
  Data *data;
};


#endif //GRADIENTFOREST_REGRESSIONSPLITTINGRULEFACTORY_H
