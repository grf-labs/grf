#ifndef RANGER_QUANTILEPREDICTIONSTRATEGY_H
#define RANGER_QUANTILEPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "PredictionStrategy.h"

class InstrumentalPredictionStrategy: public PredictionStrategy {
public:
  InstrumentalPredictionStrategy(size_t instrument_varID,
                                 size_t treatment_varID,
                                 size_t dependent_varID,
                                 std::unordered_map<size_t, std::vector<double>>* responses);
  std::vector<double> predict(std::unordered_map<size_t, double> &weights_by_sampleID);

private:
  size_t instrument_varID;
  size_t treatment_varID;
  size_t dependent_varID;
  std::unordered_map<size_t, std::vector<double>>* responses;
};


#endif //RANGER_QUANTILEPREDICTIONSTRATEGY_H
