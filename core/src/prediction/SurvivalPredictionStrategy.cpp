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

#include <cmath>
#include "prediction/SurvivalPredictionStrategy.h"

namespace grf {

SurvivalPredictionStrategy::SurvivalPredictionStrategy(size_t num_failures) :
  num_failures(num_failures) {};

size_t SurvivalPredictionStrategy::prediction_length() const {
    return num_failures;
}

std::vector<double> SurvivalPredictionStrategy::predict(const std::vector<double>& average) const {
  return average;
}

std::vector<double> SurvivalPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {
  return { 0.0 };
}

size_t SurvivalPredictionStrategy::prediction_value_length() const {
  return 1;
}

PredictionValues SurvivalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_samples.at(i);
    if (leaf_node.empty()) {
      continue;
    }

    // the failure times will always range from 0, ..., num_failures
    // where num_failures is the count of failures in the training data.
    std::vector<double> count_failure(num_failures + 1);
    std::vector<double> count_censor(num_failures + 1);
    double sum = 0;

    for (auto& sample : leaf_node) {
      size_t failure_time = data.get_outcome(sample);
      double sample_weight = data.get_weight(sample);
      if (data.get_censor(sample)) {
        count_failure[failure_time] += sample_weight;
      } else {
        count_censor[failure_time] += sample_weight;
      }
      sum += sample_weight;
    }
    // Kaplan–Meier estimator of the survival function S(t)
    double kaplan_meier = 1;
    sum = sum - count_censor[0];
    std::vector<double>& survival_function = values[i];
    survival_function.resize(num_failures);

    for (size_t time = 1; time <= num_failures; time++) {
      if (sum > 0) {
        kaplan_meier = kaplan_meier * (1 - count_failure[time] / sum);
      }
      survival_function[time - 1] = kaplan_meier;
      sum = sum - count_failure[time] - count_censor[time];
    }
  }

  return PredictionValues(values, num_failures);
}

std::vector<std::pair<double, double>> SurvivalPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
