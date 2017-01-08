#ifndef GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "Data.h"
#include "PredictionStrategy.h"

class InstrumentalPredictionStrategy: public PredictionStrategy {
public:
  std::vector<double> predict(std::unordered_map<size_t, double> &weights_by_sampleID,
                              Observations observations);
};

#endif //GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H
