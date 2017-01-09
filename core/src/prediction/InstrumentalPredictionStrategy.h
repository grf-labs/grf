#ifndef GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "Data.h"
#include "PredictionStrategy.h"

class InstrumentalPredictionStrategy: public PredictionStrategy {
public:
  size_t prediction_length();
  std::vector<double> predict(const std::unordered_map<size_t, double> &weights_by_sampleID,
                              const Observations& observations);
};

#endif //GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H
