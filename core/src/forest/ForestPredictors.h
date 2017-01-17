#ifndef GRADIENTFOREST_FORESTPREDICTORS_H
#define GRADIENTFOREST_FORESTPREDICTORS_H

#include "ForestPredictor.h"

class ForestPredictors {
public:
  static ForestPredictor instrumental_predictor(uint num_threads,
                                                uint ci_group_size);

  static ForestPredictor quantile_predictor(uint num_threads,
                                            const std::vector<double>& quantiles);

  static ForestPredictor regression_predictor(uint num_threads);
};


#endif //GRADIENTFOREST_FORESTPREDICTORS_H
