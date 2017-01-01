#include "RegressionPredictionStrategy.h"

std::vector<double> RegressionPredictionStrategy::predict(std::unordered_map<size_t, double>& weights_by_sampleID,
                                                          Observations* observations) {
  double numerator = 0.0;
  double denominator = 0.0;

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t sampleID = it->first;
    double weight = it->second;
    double outcome = observations->get(Observations::OUTCOME)[sampleID];

    numerator += outcome * weight;
    denominator += weight;
  }

  return { numerator / denominator };
}