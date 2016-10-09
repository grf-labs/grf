#include "QuantileRelabelingStrategy.h"

QuantileRelabelingStrategy::QuantileRelabelingStrategy(std::vector<double>* quantiles) :
    quantiles(quantiles) {}

std::vector<uint>* QuantileRelabelingStrategy::relabelResponses(std::vector<double>* responses) {
  std::vector<double> sorted_responses(*responses);
  std::sort(sorted_responses.begin(), sorted_responses.end());

  size_t num_samples = responses->size();
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
  std::vector<uint>* relabeled_responses = new std::vector<uint>();
  for (auto& response : *responses) {
    auto quantile = std::lower_bound(quantile_cutoffs.begin(),
                                     quantile_cutoffs.end(),
                                     response);
    long quantile_index = quantile - quantile_cutoffs.begin();
    relabeled_responses->push_back((uint) quantile_index);
  }
  return relabeled_responses;
}
