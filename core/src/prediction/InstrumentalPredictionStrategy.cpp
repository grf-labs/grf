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
#include <iostream>
#include <vector>

#include "commons/Data.h"
#include "commons/utility.h"
#include "prediction/InstrumentalPredictionStrategy.h"

const std::size_t InstrumentalPredictionStrategy::OUTCOME = 0;
const std::size_t InstrumentalPredictionStrategy::TREATMENT = 1;
const std::size_t InstrumentalPredictionStrategy::INSTRUMENT = 2;
const std::size_t InstrumentalPredictionStrategy::OUTCOME_INSTRUMENT = 3;
const std::size_t InstrumentalPredictionStrategy::TREATMENT_INSTRUMENT = 4;

const std::size_t NUM_TYPES = 5;

size_t InstrumentalPredictionStrategy::prediction_length() {
    return 1;
}

std::vector<double> InstrumentalPredictionStrategy::predict(const std::vector<double>& average) {
  double instrument_effect_numerator = average.at(OUTCOME_INSTRUMENT) - average.at(OUTCOME) * average.at(INSTRUMENT);
  double first_stage_numerator = average.at(TREATMENT_INSTRUMENT) - average.at(TREATMENT) * average.at(INSTRUMENT);

  return { instrument_effect_numerator / first_stage_numerator };
}

std::vector<double> InstrumentalPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) {

  double instrument_effect_numerator = average.at(OUTCOME_INSTRUMENT)
     - average.at(OUTCOME) * average.at(INSTRUMENT);
  double first_stage_numerator = average.at(TREATMENT_INSTRUMENT)
     - average.at(TREATMENT) * average.at(INSTRUMENT);
  double treatment_effect_estimate = instrument_effect_numerator / first_stage_numerator;
  double main_effect = average.at(OUTCOME) - average.at(TREATMENT) * treatment_effect_estimate;

  double num_good_groups = 0;
  std::vector<std::vector<double>> psi_squared = {{0, 0}, {0, 0}};
  std::vector<std::vector<double>> psi_grouped_squared = {{0, 0}, {0, 0}};

  for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_values.empty(group * ci_group_size + j)) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_psi_1 = 0;
    double group_psi_2 = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {

      size_t i = group * ci_group_size + j;
      const std::vector<double>& leaf_value = leaf_values.get_values(i);

      double psi_1 = leaf_value.at(OUTCOME_INSTRUMENT)
                     - leaf_value.at(TREATMENT_INSTRUMENT) * treatment_effect_estimate
                     - leaf_value.at(INSTRUMENT) * main_effect;
      double psi_2 = leaf_value.at(OUTCOME)
                     - leaf_value.at(TREATMENT) * treatment_effect_estimate
                     - main_effect;

      psi_squared[0][0] += psi_1 * psi_1;
      psi_squared[0][1] += psi_1 * psi_2;
      psi_squared[1][0] += psi_2 * psi_1;
      psi_squared[1][1] += psi_2 * psi_2;

      group_psi_1 += psi_1;
      group_psi_2 += psi_2;
    }

    group_psi_1 /= ci_group_size;
    group_psi_2 /= ci_group_size;

    psi_grouped_squared[0][0] += group_psi_1 * group_psi_1;
    psi_grouped_squared[0][1] += group_psi_1 * group_psi_2;
    psi_grouped_squared[1][0] += group_psi_2 * group_psi_1;
    psi_grouped_squared[1][1] += group_psi_2 * group_psi_2;
  }

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      psi_squared[i][j] /= (num_good_groups * ci_group_size);
      psi_grouped_squared[i][j] /= num_good_groups;
    }
  }

  // Using notation from the GRF paper, we want to apply equation (16),
  // \hat{sigma^2} = \xi' V^{-1} Hn V^{-1}' \xi
  // with Hn = Psi as computed above, \xi = (1 0), and
  // V(x) = (E[ZW|X=x] E[Z|X=x]; E[W|X=x] 1).
  // By simple algebra, we can verify that
  // V^{-1}'\xi = 1/(E[ZW|X=x] - E[W|X=x]E[Z|X=x]) (1; -E[Z|X=x]),
  // leading to the expression below for the variance if we
  // use forest-kernel averages to estimate all conditional moments above.

  double avg_Z = average.at(INSTRUMENT);

  double var_between = 1 / (first_stage_numerator * first_stage_numerator)
    * (psi_grouped_squared[0][0]
	     - psi_grouped_squared[0][1] * avg_Z
	     - psi_grouped_squared[1][0] * avg_Z
	     + psi_grouped_squared[1][1] * avg_Z * avg_Z);

  double var_total = 1 / (first_stage_numerator * first_stage_numerator)
    * (psi_squared[0][0]
	     - psi_squared[0][1] * avg_Z
	     - psi_squared[1][0] * avg_Z
	     + psi_squared[1][1] * avg_Z * avg_Z);

  // This is the amount by which var_between is inflated due to using small groups
  double group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an objective
  // Bayes analysis of variance instead to avoid negative values.
  double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);

  double variance_estimate = var_debiased;
  return { variance_estimate };
}

size_t InstrumentalPredictionStrategy::prediction_value_length() {
  return NUM_TYPES;
}

PredictionValues InstrumentalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data* data) {
  size_t num_leaves = leaf_samples.size();

  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < leaf_samples.size(); ++i) {
    size_t leaf_size = leaf_samples[i].size();
    if (leaf_size == 0) {
      continue;
    }

    std::vector<double>& value = values[i];
    value.resize(NUM_TYPES);

    double sum_Y = 0;
    double sum_W = 0;
    double sum_Z = 0;
    double sum_YZ = 0;
    double sum_WZ = 0;

    double sum_weight = 0.0;
    for (auto& sample : leaf_samples[i]) {
      auto weight = data->get_weight(sample);
      sum_Y +=  weight * data->get_outcome(sample);
      sum_W +=  weight * data->get_treatment(sample);
      sum_Z +=  weight * data->get_instrument(sample);
      sum_YZ += weight * data->get_outcome(sample) * data->get_instrument(sample);
      sum_WZ += weight * data->get_treatment(sample) * data->get_instrument(sample);
      sum_weight += weight;
    }

    // if total weight is very small, treat the leaf as empty
    if (std::abs(sum_weight) <= 1e-16) {
      continue;
    }

    value[OUTCOME] = sum_Y / sum_weight;
    value[TREATMENT] = sum_W / sum_weight;
    value[INSTRUMENT] = sum_Z / sum_weight;
    value[OUTCOME_INSTRUMENT] = sum_YZ / sum_weight;
    value[TREATMENT_INSTRUMENT] = sum_WZ / sum_weight;
  }
  
  return PredictionValues(values, NUM_TYPES);
}

std::vector<std::pair<double, double>> InstrumentalPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data* data) {

  double instrument_effect_numerator = average.at(OUTCOME_INSTRUMENT) - average.at(OUTCOME) * average.at(INSTRUMENT);
  double first_stage_numerator = average.at(TREATMENT_INSTRUMENT) - average.at(TREATMENT) * average.at(INSTRUMENT);
  double treatment_effect_estimate = instrument_effect_numerator / first_stage_numerator;

  double outcome = data->get_outcome(sample);
  double treatment = data->get_treatment(sample);
  double instrument = data->get_instrument(sample);

  // To justify the squared residual below as an error criterion in the case of CATE estimation
  // with an unconfounded treatment assignment, see Nie and Wager (2017).
  double residual = outcome - (treatment - average.at(TREATMENT)) * treatment_effect_estimate - average.at(OUTCOME);
  double error_raw = residual * residual;

  // Estimates the Monte Carlo bias of the raw error via the jackknife estimate of variance.
  size_t num_trees = 0;
  for (size_t n = 0; n < leaf_values.get_num_nodes(); n++) {
    if (leaf_values.empty(n)) {
      continue;
    }
    num_trees++;
  }

  // If the treatment effect estimate is due to less than 5 trees, do not attempt to estimate error,
  // as this quantity is unstable due to non-linearities.
  if (num_trees <= 5) {
    return { std::make_pair<double, double>(NAN, NAN) };
  }

  // Compute 'leave one tree out' treatment effect estimates, and use them get a jackknife estimate of the excess error.
  double error_bias = 0.0;
  for (size_t n = 0; n < leaf_values.get_num_nodes(); n++) {
    if (leaf_values.empty(n)) {
      continue;
    }
    const std::vector<double>& leaf_value = leaf_values.get_values(n);
    double outcome_loto = (num_trees *  average.at(OUTCOME) - leaf_value.at(OUTCOME)) / (num_trees - 1);
    double treatment_loto = (num_trees *  average.at(TREATMENT) - leaf_value.at(TREATMENT)) / (num_trees - 1);
    double instrument_loto = (num_trees *  average.at(INSTRUMENT) - leaf_value.at(INSTRUMENT)) / (num_trees - 1);
    double outcome_instrument_loto = (num_trees *  average.at(OUTCOME_INSTRUMENT) - leaf_value.at(OUTCOME_INSTRUMENT)) / (num_trees - 1);
    double treatment_instrument_loto = (num_trees *  average.at(TREATMENT_INSTRUMENT) - leaf_value.at(TREATMENT_INSTRUMENT)) / (num_trees - 1);
    double instrument_effect_numerator_loto = outcome_instrument_loto - outcome_loto * instrument_loto;
    double first_stage_numerator_loto = treatment_instrument_loto - treatment_loto * instrument_loto;
    double treatment_effect_estimate_loto = instrument_effect_numerator_loto / first_stage_numerator_loto;
    double residual_loto = outcome - (treatment - treatment_loto) * treatment_effect_estimate_loto - outcome_loto;
    error_bias += (residual_loto - residual) * (residual_loto - residual);
  }


  error_bias *= ((double) (num_trees - 1)) / num_trees;

  double debiased_error = error_raw - error_bias;

  auto output = std::make_pair(debiased_error, error_bias);
  return {output};

}

