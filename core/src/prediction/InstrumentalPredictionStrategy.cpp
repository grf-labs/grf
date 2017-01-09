#include <vector>
#include <string>
#include "Observations.h"
#include "utility.h"
#include "InstrumentalPredictionStrategy.h"

size_t InstrumentalPredictionStrategy::prediction_length() {
    return 1;
}

std::vector<double> InstrumentalPredictionStrategy::predict(const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                            const Observations& observations) {
  // Compute the relevant averages.e
  double average_instrument = 0.0;
  double average_treatment = 0.0;
  double average_outcome = 0.0;

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    average_outcome += weight * observations.get(Observations::OUTCOME)[neighborID];
    average_treatment += weight * observations.get(Observations::TREATMENT)[neighborID];
    average_instrument += weight * observations.get(Observations::INSTRUMENT)[neighborID];
  }

  // Finally, calculate the prediction.
  double numerator = 0.0;
  double denominator = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    double response = observations.get(Observations::OUTCOME)[neighborID];
    double treatment = observations.get(Observations::TREATMENT)[neighborID];
    double instrument = observations.get(Observations::INSTRUMENT)[neighborID];

    numerator += weight * (instrument - average_instrument) * (response - average_outcome);
    denominator += weight * (instrument - average_instrument) * (treatment - average_treatment);
  }

  return { numerator / denominator };
}