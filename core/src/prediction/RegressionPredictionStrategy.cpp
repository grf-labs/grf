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

#include <cmath>
#include "prediction/RegressionPredictionStrategy.h"

namespace grf {

const size_t RegressionPredictionStrategy::OUTCOME = 0;
const size_t RegressionPredictionStrategy::WEIGHT = 1;

size_t RegressionPredictionStrategy::prediction_length() const {
    return 1;
}

std::vector<double> RegressionPredictionStrategy::predict(const std::vector<double>& average) const {
  return { average.at(OUTCOME) / average.at(WEIGHT) };
}

/**
 * In general, the basic "bootstrap of little bags" algorithm, as described in Section 4.1
 * of the GRF paper (Athey & al, 2019) could be applied to regression forests. However,
 * when sampling weights are present, we need to use a delta-method based argument as descibed
 * in (15) and (16) of that paper. Specifically, with sample weights Gi, the estimating
 * equation psi_mu for regression and associated Hessian V are:
 *
 * E[psi_{mu(x)}(Yi, Gi)] = 0,  psi_mu(Yi, Gi) = Gi(Yi - mu),
 * V(x) = d/dmu E[psi_mu(Yi, Gi)] = E[Gi | Xi = x].
 *
 * Thus, following (16), the delta method involves applying the boostrap of little bags
 * with pseudo-outcomes
 *
 * rho_i = psi_{hat{mu}(x)}(Yi, Gi) / hat{V}(x) = Gi(Yi - hat{mu}) / hat{E[Gi | Xi = x]},
 *
 * where the required estimates are obtained via forest-weighted averaging.
 */
std::vector<double> RegressionPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {

  double average_weight = average.at(WEIGHT);
  double average_outcome = average.at(OUTCOME) / average_weight;

  double num_good_groups = 0;
  double rho_squared = 0;
  double rho_grouped_squared = 0;

  for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_values.empty(group * ci_group_size + j)) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    double group_rho = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {
      size_t i = group * ci_group_size + j;
      double rho = (leaf_values.get(i, OUTCOME) - average_outcome * leaf_values.get(i, WEIGHT)) / average_weight;

      rho_squared += rho * rho;
      group_rho += rho;
    }

    group_rho /= ci_group_size;
    rho_grouped_squared += group_rho * group_rho;
  }

  double var_between = rho_grouped_squared / num_good_groups;
  double var_total = rho_squared / (num_good_groups * ci_group_size);

  // This is the amount by which var_between is inflated due to using small groups
  double group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an objective
  // Bayes analysis of variance instead to avoid negative values.
  double var_debiased = bayes_debiaser.debias(var_between, group_noise, num_good_groups);

  return { var_debiased };
}


size_t RegressionPredictionStrategy::prediction_value_length() const {
  return 2;
}

PredictionValues RegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
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
      sum += data.get_weight(sample) * data.get_outcome(sample);
      weight += data.get_weight(sample);
    }

    // if total weight is very small, treat the leaf as empty
    if (std::abs(weight) <= 1e-16) {
      continue;
    }

    std::vector<double>& averages = values[i];
    averages.resize(2);
    averages[OUTCOME] = sum / leaf_node.size();
    averages[WEIGHT] = weight / leaf_node.size();
  }

  return PredictionValues(values, 2);
}

std::vector<std::pair<double, double>> RegressionPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  double outcome = data.get_outcome(sample);

  double average_weight = average.at(WEIGHT);
  double average_outcome = average.at(OUTCOME) / average_weight;
  double error = average_outcome - outcome;
  double mse = error * error;

  double bias = 0.0;
  size_t num_trees = 0;
  for (size_t n = 0; n < leaf_values.get_num_nodes(); n++) {
    if (leaf_values.empty(n)) {
      continue;
    }

    double tree_variance = (leaf_values.get(n, OUTCOME) - average_outcome * leaf_values.get(n, WEIGHT)) / average_weight;
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

} // namespace grf
