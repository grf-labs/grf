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

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include "Observations.h"
#include "utility.h"
#include "InstrumentalPredictionStrategy.h"

const std::size_t InstrumentalPredictionStrategy::OUTCOME = 0;
const std::size_t InstrumentalPredictionStrategy::TREATMENT = 1;
const std::size_t InstrumentalPredictionStrategy::INSTRUMENT = 2;
const std::size_t InstrumentalPredictionStrategy::OUTCOME_INSTRUMENT = 3;
const std::size_t InstrumentalPredictionStrategy::TREATMENT_INSTRUMENT = 4;

size_t InstrumentalPredictionStrategy::prediction_length() {
    return 1;
}

Prediction InstrumentalPredictionStrategy::predict(const std::vector<double>& averages,
                                                   const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                   const Observations& observations) {
  double instrument_effect = averages.at(OUTCOME_INSTRUMENT) - averages.at(OUTCOME) * averages.at(INSTRUMENT);
  double first_stage = averages.at(TREATMENT_INSTRUMENT) - averages.at(TREATMENT) * averages.at(INSTRUMENT);

  std::vector<double> prediction = { instrument_effect / first_stage };
  return Prediction(prediction);
}

Prediction InstrumentalPredictionStrategy::predict_with_variance(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations,
    uint ci_group_size) {

  size_t num_trees = leaf_sampleIDs.size();

  std::vector<double> leaf_Y(num_trees);
  std::vector<double> leaf_W(num_trees);
  std::vector<double> leaf_Z(num_trees);
  std::vector<double> leaf_YZ(num_trees);
  std::vector<double> leaf_WZ(num_trees);

  double total_Y = 0;
  double total_W = 0;
  double total_Z = 0;
  double total_YZ = 0;
  double total_WZ = 0;

  double non_empty_leaves = 0;

  for (size_t i = 0; i < leaf_sampleIDs.size(); ++i) {

    size_t leaf_size = leaf_sampleIDs[i].size();

    if (leaf_size == 0) {
      leaf_Y[i] = NAN;
      leaf_W[i] = NAN;
      leaf_Z[i] = NAN;
      leaf_YZ[i] = NAN;
      leaf_WZ[i] = NAN;
      continue;
    }

    double sum_Y = 0;
    double sum_W = 0;
    double sum_Z = 0;
    double sum_YZ = 0;
    double sum_WZ = 0;

    for (auto& sampleID : leaf_sampleIDs[i]) {
      sum_Y += observations.get(Observations::OUTCOME, sampleID);
      sum_W += observations.get(Observations::TREATMENT, sampleID);
      sum_Z += observations.get(Observations::INSTRUMENT, sampleID);
      sum_YZ += observations.get(Observations::OUTCOME, sampleID) * observations.get(Observations::INSTRUMENT, sampleID);
      sum_WZ += observations.get(Observations::TREATMENT, sampleID) * observations.get(Observations::INSTRUMENT, sampleID);
    }

    leaf_Y[i] = sum_Y / leaf_size;
    leaf_W[i] = sum_W / leaf_size;
    leaf_Z[i] = sum_Z / leaf_size;
    leaf_YZ[i] = sum_YZ / leaf_size;
    leaf_WZ[i] = sum_WZ / leaf_size;

    total_Y += sum_Y / leaf_size;
    total_W += sum_W / leaf_size;
    total_Z += sum_Z / leaf_size;
    total_YZ += sum_YZ / leaf_size;
    total_WZ += sum_WZ / leaf_size;

    non_empty_leaves++;
  }

  double avg_Y = total_Y / non_empty_leaves;
  double avg_W = total_W / non_empty_leaves;
  double avg_Z = total_Z / non_empty_leaves;
  double avg_YZ = total_YZ / non_empty_leaves;
  double avg_WZ = total_WZ / non_empty_leaves;

  double instrument_effect = avg_YZ - avg_Y * avg_Z;
  double first_stage = avg_WZ - avg_W * avg_Z;
  double treatment_estimate = instrument_effect / first_stage;
  double main_effect = avg_Y - avg_W * treatment_estimate;

  double num_good_groups = 0;
  std::vector<double> psi_mean = {0, 0}; // need to center this over good groups
  std::vector<std::vector<double>> psi_squared = {{0, 0}, {0, 0}};
  std::vector<std::vector<double>> psi_grouped_squared = {{0, 0}, {0, 0}};

  for (size_t group = 0; group < leaf_sampleIDs.size() / ci_group_size; ++group) {

    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_sampleIDs[group * ci_group_size + j].size() == 0) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_psi_1 = 0;
    double group_psi_2 = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {

      size_t i = group * ci_group_size + j;
      double psi_1 = leaf_YZ[i] - leaf_WZ[i] * treatment_estimate - leaf_Z[i] * main_effect;
      double psi_2 = leaf_Y[i] - leaf_W[i] * treatment_estimate - main_effect;

      psi_squared[0][0] += psi_1 * psi_1;
      psi_squared[0][1] += psi_1 * psi_2;
      psi_squared[1][0] += psi_2 * psi_1;
      psi_squared[1][1] += psi_2 * psi_2;

      group_psi_1 += psi_1;
      group_psi_2 += psi_2;
    }

    group_psi_1 /= ci_group_size;
    group_psi_2 /= ci_group_size;

    psi_mean[0] += group_psi_1;
    psi_mean[1] += group_psi_2;

    psi_grouped_squared[0][0] += group_psi_1 * group_psi_1;
    psi_grouped_squared[0][1] += group_psi_1 * group_psi_2;
    psi_grouped_squared[1][0] += group_psi_2 * group_psi_1;
    psi_grouped_squared[1][1] += group_psi_2 * group_psi_2;
  }

  for (size_t i = 0; i < 2; i++) {
    psi_mean[i] /= num_good_groups;
    for (size_t j = 0; j < 2; j++) {
      psi_squared[i][j] /= (num_good_groups * ci_group_size);
      psi_grouped_squared[i][j] /= num_good_groups;
    }
  }

  double var_between = psi_grouped_squared[0][0]
	  - psi_grouped_squared[0][1] * avg_W
	  - psi_grouped_squared[1][0] * avg_W
	  + psi_grouped_squared[1][1] * avg_W * avg_W;

  double var_total = psi_squared[0][0]
	  - psi_squared[0][1] * avg_W
	  - psi_squared[1][0] * avg_W
	  + psi_squared[1][1] * avg_W * avg_W;

  double within_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, but this can lead to negative variance estimates...
  double var_debiased = var_between - within_noise;
  
  // If simple variance correction is small, do an objective Bayes bias correction instead
  if (var_debiased < within_noise / 10 && num_good_groups >= 10) {

    // start by computing chi-squared log-density, taking within_noise as fixed
    double lx[100];
    double x[100];
    for (size_t iter = 0; iter < 100; iter++) {
      x[iter] = iter * within_noise / 400.0;
      lx[iter] = std::log(x[iter] + within_noise ) * (1.0 - num_good_groups / 2.0)
            - (num_good_groups * var_between) / (2 * (x[iter] + within_noise));
    }

    // compute maximal log-density, to stabilize call to exp() below
    double maxlx = lx[1];
    for (size_t iter = 0; iter < 100; iter++) {
      maxlx = std::max(maxlx, lx[iter]);
    }

    // now do Bayes rule
    double numerator = 0;
    double denominator = 0;
    for (size_t iter = 0; iter < 100; iter++) {
      double fx = std::exp(lx[iter] - maxlx);
      numerator += x[iter] * fx;
      denominator += fx;
    }
    var_debiased = numerator / denominator;
    
    // avoid jumpy behavior at the threshold
    var_debiased = std::min(var_debiased, within_noise / 10);
  }

  std::vector<double> predictions = { treatment_estimate };
  double variance_estimate = var_debiased / (first_stage * first_stage);
  std::vector<double> variance_estimates = { variance_estimate };
  return Prediction(predictions, variance_estimates);
}

bool InstrumentalPredictionStrategy::requires_leaf_sampleIDs() {
  return false;
}

size_t InstrumentalPredictionStrategy::prediction_values_length() {
  return 5;
}

PredictionValues InstrumentalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations) {
  size_t num_trees = leaf_sampleIDs.size();

  std::vector<std::vector<double>> values(prediction_values_length());
  
  values[OUTCOME].resize(num_trees);
  values[TREATMENT].resize(num_trees);
  values[INSTRUMENT].resize(num_trees);
  values[OUTCOME_INSTRUMENT].resize(num_trees);
  values[TREATMENT_INSTRUMENT].resize(num_trees);

  for (size_t i = 0; i < leaf_sampleIDs.size(); ++i) {
    size_t leaf_size = leaf_sampleIDs[i].size();

    if (leaf_size == 0) {
      values[OUTCOME][i] = NAN;
      values[TREATMENT][i] = NAN;
      values[INSTRUMENT][i] = NAN;;
      values[OUTCOME_INSTRUMENT][i] = NAN;;
      values[TREATMENT_INSTRUMENT][i] = NAN;;
      continue;
    }

    double sum_Y = 0;
    double sum_W = 0;
    double sum_Z = 0;
    double sum_YZ = 0;
    double sum_WZ = 0;

    for (auto& sampleID : leaf_sampleIDs[i]) {
      sum_Y += observations.get(Observations::OUTCOME, sampleID);
      sum_W += observations.get(Observations::TREATMENT, sampleID);
      sum_Z += observations.get(Observations::INSTRUMENT, sampleID);
      sum_YZ += observations.get(Observations::OUTCOME, sampleID) * observations.get(Observations::INSTRUMENT, sampleID);
      sum_WZ += observations.get(Observations::TREATMENT, sampleID) * observations.get(Observations::INSTRUMENT, sampleID);
    }

    values[OUTCOME][i] = sum_Y / leaf_size;
    values[TREATMENT][i] = sum_W / leaf_size;
    values[INSTRUMENT][i] = sum_Z / leaf_size;
    values[OUTCOME_INSTRUMENT][i] = sum_YZ / leaf_size;
    values[TREATMENT_INSTRUMENT][i] = sum_WZ / leaf_size;
  }
  
  return PredictionValues(values, num_trees);
}
