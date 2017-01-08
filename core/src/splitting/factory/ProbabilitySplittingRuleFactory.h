#ifndef GRADIENTFOREST_PROBABILITYSPLITTINGRULEFACTORY_H
#define GRADIENTFOREST_PROBABILITYSPLITTINGRULEFACTORY_H

#include <vector>

#include "globals.h"
#include "Data.h"
#include "SplittingRuleFactory.h"

class ProbabilitySplittingRuleFactory: public SplittingRuleFactory {
public:
  ProbabilitySplittingRuleFactory(Data* data, size_t num_classes);
  std::shared_ptr<SplittingRule> create();

private:
  Data *data;
  size_t num_classes;

  DISALLOW_COPY_AND_ASSIGN(ProbabilitySplittingRuleFactory);
};

#endif //GRADIENTFOREST_PROBABILITYSPLITTINGRULEFACTORY_H
