#ifndef GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "Observations.h"
#include "PredictionStrategy.h"
#include "PredictionValues.h"

class QuantilePredictionStrategy: public PredictionStrategy {
public:
  QuantilePredictionStrategy(std::vector<double> quantiles);

  size_t prediction_length();
  Prediction predict(const std::map<std::string, double>& average_prediction_values,
                     const std::unordered_map<size_t, double>& weights_by_sampleID,
                     const Observations& observations);

  Prediction predict_with_variance(const std::vector<std::vector<size_t>>& leaf_sampleIDs,
                                   const Observations& observations,
                                   uint ci_group_size);

  bool requires_leaf_sampleIDs();
  PredictionValues precompute_prediction_values(
      const std::vector<std::vector<size_t>>& leaf_sampleIDs,
      const Observations& observations);

private:
  std::vector<double> compute_quantile_cutoffs(const std::unordered_map<size_t, double> &weights_by_sampleID,
                                               std::vector<std::pair<size_t, double>> &sampleIDs_and_values);

  std::vector<double> quantiles;
};


#endif //GRADIENTFOREST_QUANTILEPREDICTIONSTRATEGY_H
