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
#include <numeric>
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
    const std::pair<std::vector<size_t>, std::vector<double>>& weights_by_sample,
    const Data& train_data,
    const Data& data) const {
  std::vector<double> values;
  values.reserve(weights_by_sample.first.size());
  for (size_t i = 0; i < weights_by_sample.first.size(); i++) {
    size_t sample = weights_by_sample.first[i];
    values.push_back(train_data.get_outcome(sample));
  }

  return compute_quantile_cutoffs(weights_by_sample, values);
}

std::vector<double> QuantilePredictionStrategy::compute_quantile_cutoffs(
    const std::pair<std::vector<size_t>, std::vector<double>>& weights_by_sample,
    const std::vector<double>& values) const {
  const auto& samples = weights_by_sample.first;
  const auto& weights = weights_by_sample.second;
  std::vector<size_t> sorted_index(samples.size());
  std::iota(sorted_index.begin(), sorted_index.end(), 0);
  std::sort(sorted_index.begin(),
            sorted_index.end(),
            [&](size_t first_index, size_t second_index) {
              // Note: we add a tie-breaker here to ensure that this sort consistently produces the
              // same element ordering. Otherwise, different runs of the algorithm could result in
              // different quantile predictions on the same data.
              return values[first_index] < values[second_index]
                  || (values[first_index] == values[second_index] && samples[first_index] < samples[second_index]);
            });

  std::vector<double> quantile_cutoffs;
  auto quantile_it = quantiles.begin();
  double cumulative_weight = 0.0;

  for (auto index : sorted_index) {
    size_t sample = samples[index];
    double value = values[index];

    cumulative_weight += weights[index];
    while (quantile_it != quantiles.end() && cumulative_weight >= *quantile_it) {
      quantile_cutoffs.push_back(value);
      ++quantile_it;
    }
  }

  double last_value = values[sorted_index.back()];
  for (; quantile_it != quantiles.end(); ++quantile_it) {
    quantile_cutoffs.push_back(last_value);
  }
  return quantile_cutoffs;
}

std::vector<double> QuantilePredictionStrategy::compute_variance(
    size_t sampleID,
    const std::vector<std::vector<size_t>>& samples_by_tree,
    const std::pair<std::vector<size_t>, std::vector<double>>& weights_by_sampleID,
    const Data& train_data,
    const Data& data,
    size_t ci_group_size) const {
  return { 0.0 };
}

} // namespace grf
