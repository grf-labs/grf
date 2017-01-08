#ifndef GRADIENTFOREST_REGRESSIONPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_REGRESSIONPREDICTIONSTRATEGY_H

#include "PredictionStrategy.h"

class RegressionPredictionStrategy: public PredictionStrategy {
public:
  std::vector<double> predict(std::unordered_map<size_t, double>& weights_by_sampleID,
                              Observations observations);
};


#endif //GRADIENTFOREST_REGRESSIONPREDICTIONSTRATEGY_H
