#ifndef RANGER_QUANTILEPREDICTIONSTRATEGY_H
#define RANGER_QUANTILEPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include <utility/Data.h>
#include "PredictionStrategy.h"

class InstrumentalPredictionStrategy: public PredictionStrategy {
public:
  std::vector<double> predict(std::unordered_map<size_t, double> &weights_by_sampleID,
                              std::unordered_map<std::string, std::vector<double>> original_observations);
};

#endif //RANGER_QUANTILEPREDICTIONSTRATEGY_H
