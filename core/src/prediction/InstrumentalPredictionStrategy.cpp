#include <vector>
#include <string>
#include "Observations.h"
#include "utility.h"
#include "InstrumentalPredictionStrategy.h"

size_t InstrumentalPredictionStrategy::prediction_length() {
    return 1;
}

Prediction InstrumentalPredictionStrategy::predict(const std::map<std::string, double>& average_prediction_values,
                                                   const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                   const Observations& observations) {
  // Compute the relevant averages.
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
  double instrument_effect = 0.0;
  double first_stage_effect = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t neighborID = it->first;
    double weight = it->second;

    double response = observations.get(Observations::OUTCOME)[neighborID];
    double treatment = observations.get(Observations::TREATMENT)[neighborID];
    double instrument = observations.get(Observations::INSTRUMENT)[neighborID];

    instrument_effect += weight * (instrument - average_instrument) * (response - average_outcome);
    first_stage_effect += weight * (instrument - average_instrument) * (treatment - average_treatment);
  }

  return Prediction({ instrument_effect / first_stage_effect });
}

Prediction InstrumentalPredictionStrategy::predict_with_variance(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations,
    uint ci_group_size) {

}

bool InstrumentalPredictionStrategy::requires_leaf_sampleIDs() {
  return true;
}

PredictionValues InstrumentalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>> &leaf_sampleIDs,
    const Observations &observations) {
  return PredictionValues();
}
