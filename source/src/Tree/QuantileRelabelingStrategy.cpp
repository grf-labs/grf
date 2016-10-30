#include <unordered_map>
#include "QuantileRelabelingStrategy.h"

QuantileRelabelingStrategy::QuantileRelabelingStrategy(std::vector<double>* quantiles, size_t dependent_varID) :
    quantiles(quantiles), dependent_varID(dependent_varID) {}

std::unordered_map<size_t, double> QuantileRelabelingStrategy::relabelResponses(Data* data, std::vector<size_t>& nodeSampleIDs) {

  std::vector<double> responses;
  for (auto& sampleID : nodeSampleIDs) {
    responses.push_back(data->get(sampleID, dependent_varID));
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
  std::unordered_map<size_t, double> relabeled_responses;
  for (size_t sampleID : nodeSampleIDs) {
    auto quantile = std::lower_bound(quantile_cutoffs.begin(),
                                     quantile_cutoffs.end(),
                                     data->get(sampleID, dependent_varID));
    long quantile_index = quantile - quantile_cutoffs.begin();
    relabeled_responses[sampleID] = ((uint) quantile_index);
  }
  return relabeled_responses;
}
