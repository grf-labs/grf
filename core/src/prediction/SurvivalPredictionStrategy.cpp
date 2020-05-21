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

std::vector<double> SurvivalPredictionStrategy::predict(size_t prediction_sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Data& train_data,
    const Data& data) const {
  // the event times will always range from 0, ..., num_failures
  // where num_failures is the count of failures in the training data.
  std::vector<double> count_failure(num_failures + 1);
  std::vector<double> count_censor(num_failures + 1);
  double sum = 0;
  for (const auto& entry : weights_by_sample) {
    size_t sample = entry.first;
    double forest_weight = entry.second;
    size_t failure_time = train_data.get_outcome(sample);
    double sample_weight = train_data.get_weight(sample);
    if (train_data.is_censored(sample)) {
     count_failure[failure_time] += forest_weight * sample_weight;
    } else {
     count_censor[failure_time] += forest_weight * sample_weight;
    }
    sum += forest_weight * sample_weight;
  }
  // Kaplan–Meier estimator of the survival function S(t)
  double kaplan_meier = 1;
  sum = sum - count_censor[0];
  std::vector<double> survival_function(num_failures);

  for (size_t time = 1; time <= num_failures; time++) {
   if (sum > 0) {
     kaplan_meier = kaplan_meier * (1 - count_failure[time] / sum);
     // If the estimate hits zero it will stay zero and we can break early.
     // This also prevents errors from accumulating which may yield some point estimates less than zero.
     if (kaplan_meier <= 0) {
       break;
     }
   }
   survival_function[time - 1] = kaplan_meier;
   sum = sum - count_failure[time] - count_censor[time];
  }

  return survival_function;
}

std::vector<double> SurvivalPredictionStrategy::compute_variance(
    size_t sample,
    const std::vector<std::vector<size_t>>& samples_by_tree,
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    const Data& train_data,
    const Data& data,
    size_t ci_group_size) const {
  return { 0.0 };
}

} // namespace grf
