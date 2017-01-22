#ifndef GRADIENTFOREST_FORESTTESTUTILITIES_H
#define GRADIENTFOREST_FORESTTESTUTILITIES_H

#include "ForestTrainer.h"


class ForestTestUtilities {
public:
  static void init_trainer(ForestTrainer& trainer);
  static void init_honest_trainer(ForestTrainer& trainer);
private:
  static void init_trainer(ForestTrainer& trainer, bool honesty);
};

#endif //GRADIENTFOREST_FORESTTESTUTILITIES_H
