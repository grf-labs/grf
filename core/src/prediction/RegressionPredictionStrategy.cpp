#include "RegressionPredictionStrategy.h"

size_t RegressionPredictionStrategy::prediction_length() {
    return 1;
}

std::vector<double> RegressionPredictionStrategy::predict(const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                          const Observations& observations) {
  double average = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t sampleID = it->first;
    double weight = it->second;

    double outcome = observations.get(Observations::OUTCOME)[sampleID];
    average += outcome * weight;
  }

  return { average };
}