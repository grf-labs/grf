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
#include "prediction/ProbabilityPredictionStrategy.h"

namespace grf {

ProbabilityPredictionStrategy::ProbabilityPredictionStrategy(size_t num_classes):
    num_classes(num_classes),
    num_types(num_classes + 1),
    weight_index(num_classes) {
};

size_t ProbabilityPredictionStrategy::prediction_length() const {
    return num_classes;
}

std::vector<double> ProbabilityPredictionStrategy::predict(const std::vector<double>& average) const {
  double weight_bar = average[weight_index];
  std::vector<double> predictions(num_classes);
  for (size_t cls = 0; cls < num_classes; ++cls) {
    predictions[cls] = average[cls] / weight_bar;
  }

  return predictions;
}

std::vector<double> ProbabilityPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {
  std::vector<double> variance_estimates(num_classes);
  double weight_bar = average[weight_index];
  for (size_t cls = 0; cls < num_classes; ++cls) {
    double average_outcome = average.at(cls) / weight_bar;

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
        double rho = (leaf_values.get(i, cls) - average_outcome * leaf_values.get(i, weight_index)) / weight_bar;

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
    variance_estimates[cls] = var_debiased;
  }

  return variance_estimates;
}

size_t ProbabilityPredictionStrategy::prediction_value_length() const {
  return num_types;
}

PredictionValues ProbabilityPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < num_leaves; i++) {
    const std::vector<size_t>& leaf_node = leaf_samples.at(i);
    if (leaf_node.empty()) {
      continue;
    }

    std::vector<double>& averages = values[i];
    averages.resize(num_types);
    double weight_sum = 0.0;
    for (auto& sample : leaf_node) {
      // The data Yi will be relabeled to integers {0, ..., num_classes - 1}
      size_t sample_class = data.get_outcome(sample);
      averages[sample_class] += data.get_weight(sample);
      weight_sum += data.get_weight(sample);
    }
    // if total weight is very small, treat the leaf as empty
    if (std::abs(weight_sum) <= 1e-16) {
      averages.clear();
      continue;
    }
    // store sufficient statistics in order
    // {class_counts_1, ..., class_counts_K, sum_weight}
    for (size_t cls = 0; cls < num_classes; ++cls) {
      averages[cls] = averages[cls] / leaf_node.size();
    }
    averages[weight_index] = weight_sum / leaf_node.size();
  }

  return PredictionValues(values, num_types);
}

std::vector<std::pair<double, double>> ProbabilityPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
