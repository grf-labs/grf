#include <vector>
#include <string>
#include "InstrumentalPredictionStrategy.h"

std::vector<double> InstrumentalPredictionStrategy::predict(std::unordered_map<size_t, double>& weights_by_sampleID,
                                                            std::unordered_map<std::string, std::vector<double>> original_observations) {
  // Compute the relevant averages.
  double average_instrument = 0.0;
  double average_treatment = 0.0;
  double average_response = 0.0;

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    average_instrument += weight * original_observations["instrument"][neighborID];
    average_treatment += weight * original_observations["treatment"][neighborID];
    average_response += weight * original_observations["outcome"][neighborID];
  }

  // Finally, calculate the prediction.
  double numerator = 0.0;
  double denominator = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    double instrument = original_observations["instrument"][neighborID];
    double treatment = original_observations["treatment"][neighborID];
    double response = original_observations["outcome"][neighborID];

    numerator += weight * (instrument - average_instrument) * (response - average_response);
    denominator += weight * (instrument - average_instrument) * (treatment - average_treatment);
  }

  return { numerator / denominator };
}