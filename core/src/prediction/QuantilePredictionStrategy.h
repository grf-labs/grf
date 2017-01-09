#ifndef GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "Observations.h"
#include "PredictionStrategy.h"

class QuantilePredictionStrategy: public PredictionStrategy {
public:
  QuantilePredictionStrategy(std::vector<double> quantiles);

  size_t prediction_length();
  std::vector<double> predict(const std::unordered_map<size_t, double>& weights_by_sampleID,
                              const Observations& observations);

private:
  std::vector<double> calculateQuantileCutoffs(const std::unordered_map<size_t,double>& weights_by_sampleID,
                                               std::vector<std::pair<size_t, double>>& sampleIDs_and_values);

  std::vector<double> quantiles;
};


#endif //GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H
