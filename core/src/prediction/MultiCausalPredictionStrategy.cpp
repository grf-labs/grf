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
#include "prediction/MultiCausalPredictionStrategy.h"

namespace grf {

MultiCausalPredictionStrategy::MultiCausalPredictionStrategy(size_t num_treatments) {
  this->num_treatments = num_treatments;
  this->num_types = 2 + 2 * num_treatments + num_treatments * num_treatments;
  this->weight_index = 0;
  this->Y_index = 1;
  this->W_index = 2;
  this->YW_index = 2 + num_treatments;
  this->WW_index = 2 + 2 * num_treatments;
}

size_t MultiCausalPredictionStrategy::prediction_length() const {
    return num_treatments;
}

std::vector<double> MultiCausalPredictionStrategy::predict(const std::vector<double>& average) const {
  // Re-construct the relevant data structures from the std::vector produced in `precompute_prediction_values`
  double weight_bar = average[weight_index];
  double Y_bar = average[Y_index];
  Eigen::Map<const Eigen::VectorXd> W_bar(average.data() + W_index, num_treatments);
  Eigen::Map<const Eigen::VectorXd> YW_bar(average.data() + YW_index, num_treatments);
  Eigen::Map<const Eigen::MatrixXd> WW_bar(average.data() + WW_index, num_treatments, num_treatments);

  Eigen::VectorXd theta = (WW_bar * weight_bar - W_bar * W_bar.transpose()).inverse() *
   (YW_bar * weight_bar - Y_bar * W_bar);

  return std::vector<double> (theta.data(), theta.data() + num_treatments);
}

std::vector<double> MultiCausalPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {

  double weight_bar = average[weight_index];
  double Y_bar = average[Y_index];
  Eigen::Map<const Eigen::VectorXd> W_bar(average.data() + W_index, num_treatments);
  Eigen::Map<const Eigen::VectorXd> YW_bar(average.data() + YW_index, num_treatments);
  Eigen::Map<const Eigen::MatrixXd> WW_bar(average.data() + WW_index, num_treatments, num_treatments);

  Eigen::VectorXd theta = (WW_bar * weight_bar - W_bar * W_bar.transpose()).inverse() *
   (YW_bar * weight_bar - Y_bar * W_bar);
  double main_effect = Y_bar - theta.transpose() * W_bar;

  double num_good_groups = 0;
  Eigen::MatrixXd psi_squared = Eigen::MatrixXd::Zero(num_treatments + 1, num_treatments + 1);
  Eigen::MatrixXd psi_grouped_squared = Eigen::MatrixXd::Zero(num_treatments + 1, num_treatments + 1);

  for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_values.empty(group * ci_group_size + j)) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;

    Eigen::VectorXd group_psi_1 = Eigen::VectorXd::Zero(num_treatments);
    double group_psi_2 = 0;

    for (size_t j = 0; j < ci_group_size; ++j) {

      size_t i = group * ci_group_size + j;
      const std::vector<double>& leaf_value = leaf_values.get_values(i);
      double leaf_weight = leaf_value[weight_index]; // LATER
      double leaf_Y = leaf_value[Y_index];
      Eigen::Map<const Eigen::VectorXd> leaf_W(leaf_value.data() + W_index, num_treatments);
      Eigen::Map<const Eigen::VectorXd> leaf_YW(leaf_value.data() + YW_index, num_treatments);
      Eigen::Map<const Eigen::MatrixXd> leaf_WW(leaf_value.data() + WW_index, num_treatments, num_treatments);

      Eigen::VectorXd psi_1 = leaf_YW
                             - leaf_WW * theta
                             - leaf_W * main_effect;
      double psi_2 = leaf_Y
                     - leaf_W.transpose() * theta
                     - main_effect;

      psi_squared.topLeftCorner(num_treatments, num_treatments) += psi_1 * psi_1.transpose();
      psi_squared.topRightCorner(num_treatments, 1) += psi_1 * psi_2;
      psi_squared.bottomLeftCorner(1, num_treatments) += psi_2 * psi_1.transpose();
      psi_squared.bottomRightCorner(1, 1) += Eigen::MatrixXd::Identity(1, 1) * psi_2 * psi_2;

      group_psi_1 += psi_1;
      group_psi_2 += psi_2;
    }

    group_psi_1 /= ci_group_size;
    group_psi_2 /= ci_group_size;

    psi_grouped_squared.topLeftCorner(num_treatments, num_treatments) += group_psi_1 * group_psi_1.transpose();
    psi_grouped_squared.topRightCorner(num_treatments, 1) += group_psi_1 * group_psi_2;
    psi_grouped_squared.bottomLeftCorner(1, num_treatments) += group_psi_2 * group_psi_1.transpose();
    psi_grouped_squared.bottomRightCorner(1, 1) += Eigen::MatrixXd::Identity(1, 1) * group_psi_2 * group_psi_2;
  }

  psi_squared /= (num_good_groups * ci_group_size);
  psi_grouped_squared /= num_good_groups;

  // Using notation from the GRF paper, we want to apply equation (16),
  // \hat{sigma^2} = \xi' V^{-1} Hn V^{-1}' \xi
  // with Hn = Psi as computed above, \xi selecting the diagonal, and
  // V(x) = (E[WW|X=x] E[W|X=x]; E[W|X=x] 1).
  Eigen::MatrixXd V(num_treatments + 1, num_treatments + 1);
  V.topLeftCorner(num_treatments, num_treatments) = WW_bar;
  V.topRightCorner(num_treatments, 1) = W_bar;
  V.bottomLeftCorner(1, num_treatments) = W_bar.transpose();
  V.bottomRightCorner(1, 1) = Eigen::MatrixXd::Identity(1, 1);
  V = V.inverse();
  Eigen::VectorXd var_between = (V * psi_grouped_squared * V.transpose()).diagonal();
  Eigen::VectorXd var_total = (V * psi_squared * V.transpose()).diagonal();

  // This is the amount by which var_between is inflated due to using small groups
  Eigen::VectorXd group_noise = (var_total - var_between) / (ci_group_size - 1);

  // A simple variance correction, would be to use:
  // var_debiased = var_between - group_noise.
  // However, this may be biased in small samples; we do an elementwise objective
  // Bayes analysis of variance instead to avoid negative values.
  std::vector<double> var_debiased(num_treatments);
  for (size_t i = 0; i < num_treatments; i++) {
    var_debiased[i] = bayes_debiaser.debias(var_between[i], group_noise[i], num_good_groups);
  }

  return var_debiased;
}

size_t MultiCausalPredictionStrategy::prediction_value_length() const {
  return num_types;
}

PredictionValues MultiCausalPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_samples,
    const Data& data) const {
  size_t num_leaves = leaf_samples.size();
  std::vector<std::vector<double>> values(num_leaves);

  for (size_t i = 0; i < leaf_samples.size(); ++i) {
    size_t num_samples = leaf_samples[i].size();
    if (num_samples == 0) {
      continue;
    }

    double sum_Y = 0;
    Eigen::VectorXd sum_W = Eigen::VectorXd::Zero(num_treatments);
    Eigen::VectorXd sum_YW = Eigen::VectorXd::Zero(num_treatments);
    Eigen::MatrixXd sum_WW = Eigen::MatrixXd::Zero(num_treatments, num_treatments);
    double sum_weight = 0.0;
    for (auto& sample : leaf_samples[i]) {
      double weight = data.get_weight(sample);
      double outcome = data.get_outcome(sample);
      Eigen::VectorXd treatment = data.get_treatments(sample);
      sum_Y += weight * outcome;
      sum_W += weight * treatment;
      sum_YW += weight * outcome * treatment;
      sum_WW += weight * treatment * treatment.transpose();
      sum_weight += weight;
    }
    // if total weight is very small, treat the leaf as empty
    if (std::abs(sum_weight) <= 1e-16) {
      continue;
    }
    std::vector<double>& value = values[i];
    value.reserve(num_types);
    // store sufficient statistics in order
    // {sum_weight, sum_Y, sum_W, sum_YW, sum_WW}

    value.push_back(sum_weight / num_samples);
    value.push_back(sum_Y / num_samples);
    for (size_t j = 0; j < num_treatments; j++) {
      value.push_back(sum_W[j] / num_samples);
    }
    for (size_t j = 0; j < num_treatments; j++) {
      value.push_back(sum_YW[j] / num_samples);
    }
    for (size_t j = 0; j < num_treatments * num_treatments; j++) {
      value.push_back(sum_WW.data()[j] / num_samples);
    }

  }

  return PredictionValues(values, num_types);
}

std::vector<std::pair<double, double>> MultiCausalPredictionStrategy::compute_error(
    size_t sample,
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    const Data& data) const {
  return { std::make_pair<double, double>(NAN, NAN) };
}

} // namespace grf
