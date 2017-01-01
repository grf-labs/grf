#include "NoopRelabelingStrategy.h"

std::unordered_map<size_t, double> NoopRelabelingStrategy::relabel_outcomes(
    Observations *observations,
    std::vector<size_t> &node_sampleIDs) {

  std::unordered_map<size_t, double> relabeled_observations;
  for (size_t sampleID : node_sampleIDs) {
    double outcome = observations->get(Observations::OUTCOME)[sampleID];
    relabeled_observations[sampleID] = outcome;
  }
  return relabeled_observations;
}
