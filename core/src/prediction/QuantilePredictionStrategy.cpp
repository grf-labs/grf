#include <vector>
#include <string>
#include "Data.h"
#include "QuantilePredictionStrategy.h"

QuantilePredictionStrategy::QuantilePredictionStrategy(std::vector<double>* quantiles):
    quantiles(quantiles) {
};

std::vector<double> QuantilePredictionStrategy::predict(std::unordered_map<size_t, double>& weights_by_sampleID,
                                                        std::unordered_map<std::string, std::vector<double>> original_observations) {
  std::vector<std::pair<size_t, double>> sampleIDs_and_values;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t sampleID = it->first;
    sampleIDs_and_values.push_back(std::pair<size_t, double>(
        sampleID, original_observations["outcome"][sampleID]));
  }

  return calculateQuantileCutoffs(weights_by_sampleID, sampleIDs_and_values);
}

std::vector<double> QuantilePredictionStrategy::calculateQuantileCutoffs(
    std::unordered_map<size_t,double> &weights_by_sampleID,
    std::vector<std::pair<size_t, double>> sampleIDs_and_values) {
  std::sort(sampleIDs_and_values.begin(),
            sampleIDs_and_values.end(),
            [](std::pair<size_t, double> first_pair, std::pair<size_t, double> second_pair) {
              return first_pair.second < second_pair.second;
            });

  std::vector<double> quantile_cutoffs;
  auto quantile_it = quantiles->begin();
  double cumulative_weight = 0.0;

  for (auto it = sampleIDs_and_values.begin(); it != sampleIDs_and_values.end()
                                               && quantile_it != quantiles->end(); ++it) {
    size_t sampleID = it->first;
    double value = it->second;

    cumulative_weight += weights_by_sampleID[sampleID];
    if (cumulative_weight >= *quantile_it) {
      quantile_cutoffs.push_back(value);
      ++quantile_it;
    }
  }

  double last_value = sampleIDs_and_values.back().second;
  for (; quantile_it != quantiles->end(); ++quantile_it) {
    quantile_cutoffs.push_back(last_value);
  }
  return quantile_cutoffs;
}