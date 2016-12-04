#ifndef RANGER_QUANTILEPREDICTIONSTRATEGY_H
#define RANGER_QUANTILEPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "PredictionStrategy.h"

class QuantilePredictionStrategy: public PredictionStrategy {
public:
  QuantilePredictionStrategy(std::vector<double>* quantiles, std::vector<double>* original_responses);
  std::vector<double> predict(std::unordered_map<size_t, double> &weights_by_sampleID);
  std::vector<double> calculateQuantileCutoffs(std::unordered_map<size_t,double> &weights_by_sampleID,
                                               std::vector<std::pair<size_t, double>> sampleIDs_and_values);

private:
  std::vector<double>* quantiles;
  std::vector<double>* original_responses;
};


#endif //RANGER_QUANTILEPREDICTIONSTRATEGY_H
