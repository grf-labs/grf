/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <algorithm>
#include <vector>
#include <string>

#include "commons/DefaultData.h"
#include "prediction/QuantilePredictionStrategy.h"

QuantilePredictionStrategy::QuantilePredictionStrategy(std::vector<double> quantiles):
    quantiles(quantiles) {
};

size_t QuantilePredictionStrategy::prediction_length() {
    return quantiles.size();
}

std::vector<double> QuantilePredictionStrategy::predict(
    size_t prediction_sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Observations& observations) {
  std::vector<std::pair<size_t, double>> samples_and_values;
  for (auto it = weights_by_sample.begin(); it != weights_by_sample.end(); it++) {
    size_t sample = it->first;
    samples_and_values.push_back(std::pair<size_t, double>(
        sample, observations.get(Observations::OUTCOME, sample)));
  }

  return compute_quantile_cutoffs(weights_by_sample, samples_and_values);
}

std::vector<double> QuantilePredictionStrategy::compute_quantile_cutoffs(
    const std::unordered_map<size_t, double>& weights_by_sample,
    std::vector<std::pair<size_t, double>>& samples_and_values) {
  std::sort(samples_and_values.begin(),
            samples_and_values.end(),
            [](std::pair<size_t, double> first_pair, std::pair<size_t, double> second_pair) {
              return first_pair.second < second_pair.second;
            });

  std::vector<double> quantile_cutoffs;
  auto quantile_it = quantiles.begin();
  double cumulative_weight = 0.0;

  for (auto it = samples_and_values.begin(); it != samples_and_values.end(); ++it) {
    size_t sample = it->first;
    double value = it->second;

    cumulative_weight += weights_by_sample.at(sample);
    while (quantile_it != quantiles.end() && cumulative_weight >= *quantile_it) {
      quantile_cutoffs.push_back(value);
      ++quantile_it;
    }
  }

  double last_value = samples_and_values.back().second;
  for (; quantile_it != quantiles.end(); ++quantile_it) {
    quantile_cutoffs.push_back(last_value);
  }
  return quantile_cutoffs;
}
