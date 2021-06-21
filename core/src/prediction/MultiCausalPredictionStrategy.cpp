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

MultiCausalPredictionStrategy::MultiCausalPredictionStrategy(size_t num_treatments,
                                                             size_t num_outcomes) {
  this->num_treatments = num_treatments;
  this->num_outcomes = num_outcomes;
  this->num_types = num_treatments * (num_treatments + num_outcomes + 1) + num_outcomes + 1;
  this->weight_index = 0;
  this->Y_index = 1;
  this->W_index = this->Y_index + num_outcomes;
  this->YW_index = this->W_index + num_treatments;
  this->WW_index = this->YW_index + num_treatments * num_outcomes;
}

size_t MultiCausalPredictionStrategy::prediction_length() const {
    return num_treatments * num_outcomes;
}

std::vector<double> MultiCausalPredictionStrategy::predict(const std::vector<double>& average) const {
  // Re-construct the relevant data structures from the std::vector produced in `precompute_prediction_values`
  double weight_bar = average[weight_index];
  Eigen::Map<const Eigen::VectorXd> Y_bar(average.data() + Y_index, num_outcomes);
  Eigen::Map<const Eigen::VectorXd> W_bar(average.data() + W_index, num_treatments);
  Eigen::Map<const Eigen::MatrixXd> YW_bar(average.data() + YW_index, num_treatments, num_outcomes);
  Eigen::Map<const Eigen::MatrixXd> WW_bar(average.data() + WW_index, num_treatments, num_treatments);

  // Why we do not do a `...ldlt().solve(...)` instead of inverse: because dim(W) is low (intended use-case)
  // and ^-1 replicates the behavior of InstrumentalPredictionStrategy for dim(W) = 1.
  Eigen::MatrixXd theta = (WW_bar * weight_bar - W_bar * W_bar.transpose()).inverse() *
   (YW_bar * weight_bar - W_bar * Y_bar.transpose());

  return std::vector<double> (theta.data(), theta.data() + prediction_length());
}

/**
 * The following calculations follow directly from {@link InstrumentalPredictionStrategy}.
 * Yi is real valued, Wi is a 1xK vector, Gi are sample weights.
 * The scoring function is:
 * psi_{mu, tau}^1(Yi, Wi, Gi) = Gi Wi (Yi - Wi tau - mu),
 * psi_{mu, tau}^2(Yi, Wi, Gi) = Gi (Yi - Wi tau - mu).
 *
 * The Hessian V(x) is:
 * V11(x) = A; V12(x) = b
 * V21(x) = b'; V22(x) = c
 * where
 * A = E[Gi Wi Wi' | Xi = x]
 * b = E[Gi Wi | Xi = x]
 * c = E[Gi | Xi = x]
 *
 * Using the matrix inverse formula in block form we get
 * rhoi = \xi' hat{V}^{-1} psi_{hat{mu}, hat{tau}}(Yi, Wi, Gi)
 *      = term1 * psi_1 - term2 * psi_2
 * where \xi selects the first K-subvector and
 * term1 = A^{-1} + 1/k A^{-1} b b' A^{-1},
 * term2 = 1/k A^{-1} b,
 * and k = c - b' A^{-1} b
 */
std::vector<double> MultiCausalPredictionStrategy::compute_variance(
    const std::vector<double>& average,
    const PredictionValues& leaf_values,
    size_t ci_group_size) const {
  if (num_outcomes > 1) {
    throw std::runtime_error("Pointwise variance estimates are only implemented for one outcome.");
  }

  double weight_bar = average[weight_index];
  double Y_bar = average[Y_index];
  Eigen::Map<const Eigen::VectorXd> W_bar(average.data() + W_index, num_treatments);
  Eigen::Map<const Eigen::VectorXd> YW_bar(average.data() + YW_index, num_treatments);
  Eigen::Map<const Eigen::MatrixXd> WW_bar(average.data() + WW_index, num_treatments, num_treatments);

  Eigen::VectorXd theta = (WW_bar * weight_bar - W_bar * W_bar.transpose()).inverse() *
   (YW_bar * weight_bar - Y_bar * W_bar);
  double main_effect = (Y_bar - theta.transpose() * W_bar) / weight_bar;

  // NOTE: could potentially use pseudoinverses.
  double k = weight_bar - W_bar.transpose() * WW_bar.inverse() * W_bar;
  Eigen::MatrixXd term1 = WW_bar.inverse() + 1 / k * WW_bar.inverse() * W_bar * W_bar.transpose() * WW_bar.inverse();
  Eigen::VectorXd term2 = 1 / k * WW_bar.inverse() * W_bar;

  double num_good_groups = 0;
  Eigen::VectorXd rho_squared = Eigen::VectorXd::Zero(num_treatments);
  Eigen::VectorXd rho_grouped_squared = Eigen::VectorXd::Zero(num_treatments);

  Eigen::VectorXd group_rho = Eigen::VectorXd(num_treatments);
  Eigen::VectorXd psi_1 = Eigen::VectorXd(num_treatments);
  Eigen::VectorXd rho = Eigen::VectorXd(num_treatments);
  for (size_t group = 0; group < leaf_values.get_num_nodes() / ci_group_size; ++group) {
    bool good_group = true;
    for (size_t j = 0; j < ci_group_size; ++j) {
      if (leaf_values.empty(group * ci_group_size + j)) {
        good_group = false;
      }
    }
    if (!good_group) continue;

    num_good_groups++;
    group_rho.setZero();
    for (size_t j = 0; j < ci_group_size; ++j) {

      size_t i = group * ci_group_size + j;
      const std::vector<double>& leaf_value = leaf_values.get_values(i);
      double leaf_weight = leaf_value[weight_index];
      double leaf_Y = leaf_value[Y_index];
      Eigen::Map<const Eigen::VectorXd> leaf_W(leaf_value.data() + W_index, num_treatments);
      Eigen::Map<const Eigen::VectorXd> leaf_YW(leaf_value.data() + YW_index, num_treatments);
      Eigen::Map<const Eigen::MatrixXd> leaf_WW(leaf_value.data() + WW_index, num_treatments, num_treatments);

      psi_1 = leaf_YW - leaf_WW * theta - leaf_W * main_effect;
      double psi_2 = leaf_Y - leaf_W.transpose() * theta - leaf_weight * main_effect;

      rho = term1 * psi_1 - term2 * psi_2;
      rho_squared += rho.array().square().matrix();
      group_rho += rho;
    }

    group_rho /= ci_group_size;
    rho_grouped_squared += group_rho.array().square().matrix();
  }

  Eigen::VectorXd var_between = rho_grouped_squared / num_good_groups;
  Eigen::VectorXd var_total = rho_squared / (num_good_groups * ci_group_size);

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

    Eigen::VectorXd sum_Y = Eigen::VectorXd::Zero(num_outcomes);
    Eigen::VectorXd sum_W = Eigen::VectorXd::Zero(num_treatments);
    Eigen::MatrixXd sum_YW = Eigen::MatrixXd::Zero(num_treatments, num_outcomes);
    Eigen::MatrixXd sum_WW = Eigen::MatrixXd::Zero(num_treatments, num_treatments);
    double sum_weight = 0.0;
    for (auto& sample : leaf_samples[i]) {
      double weight = data.get_weight(sample);
      Eigen::VectorXd outcome = data.get_outcomes(sample);
      Eigen::VectorXd treatment = data.get_treatments(sample);
      sum_Y += weight * outcome;
      sum_W += weight * treatment;
      sum_YW.noalias() += weight * treatment * outcome.transpose();
      sum_WW.noalias() += weight * treatment * treatment.transpose();
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
    for (size_t j = 0; j < num_outcomes; j++) {
      value.push_back(sum_Y[j] / num_samples);
    }
    for (size_t j = 0; j < num_treatments; j++) {
      value.push_back(sum_W[j] / num_samples);
    }
    for (size_t j = 0; j < num_treatments * num_outcomes; j++) {
      value.push_back(sum_YW.data()[j] / num_samples);
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
