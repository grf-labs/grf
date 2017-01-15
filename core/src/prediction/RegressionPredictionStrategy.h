#ifndef GRADIENTFOREST_REGRESSIONPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_REGRESSIONPREDICTIONSTRATEGY_H

#include "Data.h"
#include "Observations.h"
#include "PredictionStrategy.h"
#include "PredictionValues.h"

class RegressionPredictionStrategy: public PredictionStrategy {
public:
  size_t prediction_length();

  std::vector<double> predict(const std::map<std::string, double>& average_prediction_values,
                              const std::unordered_map<size_t, double>& weights_by_sampleID,
                              const Observations& observations);

  bool requires_leaf_sampleIDs();
  PredictionValues precompute_prediction_values(const std::vector<std::vector<size_t>>& leaf_sampleIDs,
                                                const Observations& observations);
};


#endif //GRADIENTFOREST_REGRESSIONPREDICTIONSTRATEGY_H
