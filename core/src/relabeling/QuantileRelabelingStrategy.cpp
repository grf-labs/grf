
#include <unordered_map>
#include "QuantileRelabelingStrategy.h"

QuantileRelabelingStrategy::QuantileRelabelingStrategy(std::vector<double>* quantiles) :
    quantiles(quantiles) {}

std::unordered_map<size_t, double> QuantileRelabelingStrategy::relabel_outcomes(
    Observations *observations,
    std::vector<size_t> &node_sampleIDs) {

  std::vector<double> responses;
  for (auto& sampleID : node_sampleIDs) {
    responses.push_back(observations->get(Observations::OUTCOME)[sampleID]);
  }

  std::vector<double> sorted_responses(responses);
  std::sort(sorted_responses.begin(), sorted_responses.end());

  size_t num_samples = responses.size();
  std::vector<double> quantile_cutoffs;

  // Calculate the response value cutoffs for each quantile.
  for (auto& quantile : *quantiles) {
    size_t response_index = (size_t) ceil(num_samples * quantile) - 1;
    quantile_cutoffs.push_back(sorted_responses[response_index]);
  }

  // Remove duplicate cutoffs.
  quantile_cutoffs.erase(std::unique(quantile_cutoffs.begin(), quantile_cutoffs.end()),
                         quantile_cutoffs.end());

  // Assign a class to each response based on what quantile it belongs to.
  std::unordered_map<size_t, double> relabeled_observations;
  for (size_t sampleID : node_sampleIDs) {
    double outcome = observations->get(Observations::OUTCOME)[sampleID];
    auto quantile = std::lower_bound(quantile_cutoffs.begin(),
                                     quantile_cutoffs.end(),
                                     outcome);
    long quantile_index = quantile - quantile_cutoffs.begin();
    relabeled_observations[sampleID] = ((uint) quantile_index);
  }
  return relabeled_observations;
}
