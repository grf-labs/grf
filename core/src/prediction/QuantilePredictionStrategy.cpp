/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <algorithm>
#include <vector>
#include <string>

#include "commons/Data.h"
#include "prediction/QuantilePredictionStrategy.h"

QuantilePredictionStrategy::QuantilePredictionStrategy(std::vector<double> quantiles):
    quantiles(quantiles) {
};

size_t QuantilePredictionStrategy::prediction_length() {
    return quantiles.size();
}

Prediction QuantilePredictionStrategy::predict(size_t sampleID,
                                               const std::unordered_map<size_t, double>& weights_by_sampleID,
                                               const Observations& observations) {
  std::vector<std::pair<size_t, double>> sampleIDs_and_values;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    size_t sampleID = it->first;
    sampleIDs_and_values.push_back(std::pair<size_t, double>(
        sampleID, observations.get(Observations::OUTCOME, sampleID)));
  }

  std::vector<double> quantile_cutoffs = compute_quantile_cutoffs(weights_by_sampleID, sampleIDs_and_values);
  return Prediction(quantile_cutoffs);
}

std::vector<double> QuantilePredictionStrategy::compute_quantile_cutoffs(
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    std::vector<std::pair<size_t, double>>& sampleIDs_and_values) {
  std::sort(sampleIDs_and_values.begin(),
            sampleIDs_and_values.end(),
            [](std::pair<size_t, double> first_pair, std::pair<size_t, double> second_pair) {
              return first_pair.second < second_pair.second;
            });

  std::vector<double> quantile_cutoffs;
  auto quantile_it = quantiles.begin();
  double cumulative_weight = 0.0;

  for (auto it = sampleIDs_and_values.begin(); it != sampleIDs_and_values.end(); ++it) {
    size_t sampleID = it->first;
    double value = it->second;

    cumulative_weight += weights_by_sampleID.at(sampleID);
    while (quantile_it != quantiles.end() && cumulative_weight >= *quantile_it) {
      quantile_cutoffs.push_back(value);
      ++quantile_it;
    }
  }

  double last_value = sampleIDs_and_values.back().second;
  for (; quantile_it != quantiles.end(); ++quantile_it) {
    quantile_cutoffs.push_back(last_value);
  }
  return quantile_cutoffs;
}
