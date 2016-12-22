#ifndef GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "Observations.h"
#include "PredictionStrategy.h"

class QuantilePredictionStrategy: public PredictionStrategy {
public:
  QuantilePredictionStrategy(std::vector<double>* quantiles);
  std::vector<double> predict(std::unordered_map<size_t, double>& weights_by_sampleID,
                              Observations observations);

private:
  std::vector<double> calculateQuantileCutoffs(std::unordered_map<size_t,double> &weights_by_sampleID,
                                               std::vector<std::pair<size_t, double>> sampleIDs_and_values);

  std::vector<double>* quantiles;
};


#endif //GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H
