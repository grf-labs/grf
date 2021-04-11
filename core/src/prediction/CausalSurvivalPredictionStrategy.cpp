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
#include <vector>

#include "commons/Data.h"
#include "commons/utility.h"
#include "prediction/CausalSurvivalPredictionStrategy.h"

namespace grf {

const std::size_t CausalSurvivalPredictionStrategy::NUMERATOR = 0;
const std::size_t CausalSurvivalPredictionStrategy::DENOMINATOR = 1;

const std::size_t NUM_TYPES = 2;

size_t CausalSurvivalPredictionStrategy::prediction_length() const {
  return 1;
}

std::vector<double> CausalSurvivalPredictionStrategy::predict(const std::vector<double>& average) const {
  return { average.at(NUMERATOR) / average.at(DENOMINATOR) };
}

std::vector<double> CausalSurvivalPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {

  double v_est = average.at(DENOMINATOR);
  double average_eta = average.at(NUMERATOR) / average.at(DENOMINATOR);


  double num_good_groups = 0;
  double psi_squared = 0;
  double psi_grouped_squared = 0;

  for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_values.empty(group * ci_group_size + j)) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_psi = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {

      size_t i = group * ci_group_size + j;
      const std::vector<double>& leaf_value = leaf_values.get_values(i);

      double psi_1 = leaf_value.at(NUMERATOR) - leaf_value.at(DENOMINATOR) * average_eta;

      psi_squared += psi_1 * psi_1;
      group_psi += psi_1;
    }

    group_psi /= ci_group_size;
		psi_grouped_squared += group_psi * group_psi;
  }

  // Using notation from the GRF paper, ...

  double var_between = psi_grouped_squared / num_good_groups;
  double var_total = psi_squared / (num_good_groups * ci_group_size);

  // This is the amount by which var_between is inflated due to using small groups
  double group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an objective
  // Bayes analysis of variance instead to avoid negative values.
  double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);

  double variance_estimate = var_debiased / (v_est * v_est);
  return { variance_estimate };
}

size_t CausalSurvivalPredictionStrategy::prediction_value_length() const {
  return NUM_TYPES;
}

PredictionValues CausalSurvivalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();

  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < leaf_samples.size(); ++i) {
    size_t leaf_size = leaf_samples[i].size();
    if (leaf_size == 0) {
      continue;
    }

    double numerator_sum = 0;
    double denominator_sum = 0;
    double sum_weight = 0;

    for (auto& sample : leaf_samples[i]) {
      double weight = data.get_weight(sample);
      numerator_sum += weight * data.get_causal_survival_numerator(sample);
      denominator_sum += weight * data.get_causal_survival_denominator(sample);
      sum_weight += weight;
    }

    // if total weight is very small, treat the leaf as empty
    if (std::abs(sum_weight) <= 1e-16) {
      continue;
    }
    std::vector<double>& value = values[i];
    value.resize(NUM_TYPES);

    value[NUMERATOR] = numerator_sum / leaf_size;
    value[DENOMINATOR] = denominator_sum / leaf_size;
  }

  return PredictionValues(values, NUM_TYPES);
}

std::vector<std::pair<double, double>> CausalSurvivalPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
