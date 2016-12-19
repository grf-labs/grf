#include "InstrumentalRelabelingStrategy.h"

std::unordered_map<size_t, double> InstrumentalRelabelingStrategy::relabelObservations(
    std::unordered_map<std::string, std::vector<double>> *observations,
    std::vector<size_t> &node_sampleIDs) {

  // Calculate the relevant averages;
  size_t num_samples = node_sampleIDs.size();

  double total_treatment = 0.0;
  double total_instrument = 0.0;
  double total_response = 0.0;

  for (size_t sampleID : node_sampleIDs) {
    total_treatment += (*observations)["treatment"][sampleID];
    total_instrument += (*observations)["instrument"][sampleID];
    total_response += (*observations)["outcome"][sampleID];
  }

  double average_treatment = total_treatment / num_samples;
  double average_instrument = total_instrument / num_samples;
  double average_response = total_response / num_samples;

  // Calculate the treatment effect.
  double numerator = 0.0;
  double denominator = 0.0;
  double regularizer = 0.0;
  for (size_t sampleID : node_sampleIDs) {
    double treatment = (*observations)["treatment"][sampleID];
    double instrument = (*observations)["instrument"][sampleID];
    double response = (*observations)["outcome"][sampleID];

    numerator += (instrument - average_instrument) * (response - average_response);
    denominator += (instrument - average_instrument) * (treatment - average_treatment);
    regularizer += (instrument - average_instrument) * (instrument - average_instrument);
  }

  if (equalDoubles(denominator, 0.0)) {
    return std::unordered_map<size_t, double>(); // Signals that we should not perform a split.
  }

  // 50 is a hack
  double local_average_treatment_effect = numerator / (denominator + 50 / node_sampleIDs.size() * regularizer * sgn(denominator));

  // Create the new responses;
  std::unordered_map<size_t, double> new_responses_by_sampleID;

  for (size_t sampleID : node_sampleIDs) {
    double treatment = (*observations)["treatment"][sampleID];
    double instrument = (*observations)["instrument"][sampleID];
    double response = (*observations)["outcome"][sampleID];

    double residual = (response - average_response) - local_average_treatment_effect * (treatment - average_treatment);
    new_responses_by_sampleID[sampleID] = (instrument - average_instrument) * residual;
  }
  return new_responses_by_sampleID;
}

bool InstrumentalRelabelingStrategy::equalDoubles(double first, double second) {
  return std::abs(first - second) < std::numeric_limits<double>::epsilon();
}

int InstrumentalRelabelingStrategy::sgn(double val) {
  return (0 < val) - (val < 0);
}