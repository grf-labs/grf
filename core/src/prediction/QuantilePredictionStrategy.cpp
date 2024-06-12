/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#include "commons/Data.h"
#include "prediction/QuantilePredictionStrategy.h"

namespace grf {

QuantilePredictionStrategy::QuantilePredictionStrategy(std::vector<double> quantiles):
    quantiles(quantiles) {
};

size_t QuantilePredictionStrategy::prediction_length() const {
    return quantiles.size();
}

std::vector<double> QuantilePredictionStrategy::predict(
    size_t prediction_sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Data& train_data,
    const Data& data) const {
  std::vector<std::pair<size_t, double>> samples_and_values;
  for (const auto& entry : weights_by_sample) {
    size_t sample = entry.first;
    samples_and_values.emplace_back(sample, train_data.get_outcome(sample));
  }

  return compute_quantile_cutoffs(weights_by_sample, samples_and_values);
}

std::vector<double> QuantilePredictionStrategy::compute_quantile_cutoffs(
    const std::unordered_map<size_t, double>& weights_by_sample,
    std::vector<std::pair<size_t, double>>& samples_and_values) const {
  std::sort(samples_and_values.begin(),
            samples_and_values.end(),
            [](std::pair<size_t, double> first_pair, std::pair<size_t, double> second_pair) {
              // Note: we add a tie-breaker here to ensure that this sort consistently produces the
              // same element ordering. Otherwise, different runs of the algorithm could result in
              // different quantile predictions on the same data.
              return first_pair.second < second_pair.second
                  || (first_pair.second == second_pair.second && first_pair.first < second_pair.first);
            });

  std::vector<double> quantile_cutoffs;
  auto quantile_it = quantiles.begin();
  double cumulative_weight = 0.0;

  for (const auto& entry : samples_and_values) {
    size_t sample = entry.first;
    double value = entry.second;

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

std::vector<double> QuantilePredictionStrategy::compute_variance(
    size_t sampleID,
    const std::vector<std::vector<size_t>>& samples_by_tree,
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    const Data& train_data,
    const Data& data,
    size_t ci_group_size) const {
  return { 0.0 };
}

} // namespace grf
