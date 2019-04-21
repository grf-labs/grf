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
#include <string>
#include "prediction/RegressionPredictionStrategy.h"

const size_t RegressionPredictionStrategy::OUTCOME = 0;

size_t RegressionPredictionStrategy::prediction_length() {
    return 1;
}

std::vector<double> RegressionPredictionStrategy::predict(const std::vector<double>& average) {
  return { average.at(OUTCOME) };
}

std::vector<double> RegressionPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) {

  double average_outcome = average.at(OUTCOME);

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
      double psi_1 = leaf_values.get(i, OUTCOME) - average_outcome;

      psi_squared += psi_1 * psi_1;
      group_psi += psi_1;
    }

    group_psi /= ci_group_size;
    psi_grouped_squared += group_psi * group_psi;
  }

  double var_between = psi_grouped_squared / num_good_groups;
  double var_total = psi_squared / (num_good_groups * ci_group_size);

  // This is the amount by which var_between is inflated due to using small groups
  double group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an objective
  // Bayes analysis of variance instead to avoid negative values.
  double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);

  return { var_debiased };
}


size_t RegressionPredictionStrategy::prediction_value_length() {
  return 1;
}

PredictionValues RegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data* data) {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_samples.at(i);
    if (leaf_node.empty()) {
      continue;
    }

    double sum = 0.0;
    double weight = 0.0;
    for (auto& sample : leaf_node) {
      sum += data->get_weight(sample) * data->get_outcome(sample);
      weight  += data->get_weight(sample);
    }

    // if total weight is very small, treat the leaf as empty
    if (std::abs(weight) <= 1e-16) {
      continue;
    }

    std::vector<double> &averages = values[i];
    averages.resize(1);
    averages[OUTCOME] = sum / weight;
  }

  return PredictionValues(values, 1);
}

std::vector<std::pair<double, double>>  RegressionPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data* data) {
  double outcome = data->get_outcome(sample);

  double error = average.at(OUTCOME) - outcome;
  double mse = error * error;

  double bias = 0.0;
  size_t num_trees = 0;
  for (size_t n = 0; n < leaf_values.get_num_nodes(); n++) {
    if (leaf_values.empty(n)) {
      continue;
    }

    double tree_variance = leaf_values.get(n, OUTCOME) - average.at(OUTCOME);
    bias += tree_variance * tree_variance;
    num_trees++;
  }

  if (num_trees <= 1) {
    return { std::make_pair<double, double>(NAN, NAN) };
  }

  bias /= num_trees * (num_trees - 1);

  double debiased_error = mse - bias;

  auto output = std::pair<double, double>(debiased_error, bias);
  return { output };
}

