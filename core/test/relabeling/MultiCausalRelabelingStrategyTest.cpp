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

#include "catch.hpp"
#include "commons/utility.h"
#include "relabeling/RelabelingStrategy.h"
#include "relabeling/MultiCausalRelabelingStrategy.h"

using namespace grf;

Eigen::ArrayXXd get_relabeled_outcomes(
  std::vector<double> observations, size_t num_samples, size_t num_treatments, const std::vector<double>& gradient_weights) {

  Data data(observations, num_samples, 1 + num_treatments);
  data.set_outcome_index(0);
  std::vector<size_t> treatment_index(num_treatments);
  std::iota(treatment_index.begin(), treatment_index.end(), 1);
  data.set_treatment_index(treatment_index);

  std::vector<size_t> samples;
  for (size_t i = 0; i < num_samples; ++i) {
    samples.push_back(i);
  }

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new MultiCausalRelabelingStrategy(num_treatments, gradient_weights));

  Eigen::ArrayXXd relabeled_observations(num_samples, num_treatments);
  bool stop = relabeling_strategy->relabel(samples, data, relabeled_observations);
  if (stop) {
    return Eigen::ArrayXXd(0, 0);
  }

  return relabeled_observations;
}

TEST_CASE("multi causal relabeling calculations are correct", "[multi causal, relabeling]") {
  /* This test case data is generated from the following R script
    set.seed(123)
    n <- 10
    W <- round(matrix(rnorm(n * 2), n, 2), 2)
    Y <- round(W[, 1] + 2 * W[, 2] + rnorm(n), 2)
    Y <- Y - mean(Y)
    W <- sweep(W, 2, colMeans(W), "-")

    WWinv <- solve(t(W) %*% W)
    b <- WWinv %*% t(W) %*% Y
    residuals <- Y - W %*% b
    influence <- matrix(NA, 2, n)
    for (i in 1:n) {
      influence[, i] <- WWinv %*% W[i, ] %*% residuals[i]
    }
    # `std::vector<double> observations`:
    dput(c(Y, W))
    # `Eigen::ArrayXXd rho_expected`:
    dput(influence)
  */
  size_t num_samples = 10;
  size_t num_treatments = 2;
  std::vector<double> observations = {
    0.747, 0.207, 1.267, -0.503, -1.683, 3.547, 2.237, -5.123,
    -0.493, -0.203, -0.634, -0.304, 1.486, -0.00399999999999999,
    0.056, 1.646, 0.386, -1.344, -0.764, -0.524, 1.012, 0.152, 0.192,
    -0.098, -0.768, 1.582, 0.292, -2.178, 0.492, -0.678
  }; // [Y, W_1, W_2]

  Eigen::ArrayXXd rho_expected(num_samples, num_treatments);
  rho_expected <<
    0.0573262020170922, -0.058379510635618, -0.0168308217417549,
    0.0127092292843973, -0.0871238862874838, 0.0396239647273461,
    -0.00305980821703425, 0.00513039856113015, -0.045747755401616,
    0.0668820374735852, -0.0540282235976447, -0.0311225877163698,
    0.0607359331696434, 0.0109600361706618, 0.018357949680282, 0.112964940129093,
    0.109475829631935, -0.0878823962298169, -0.0391054192534192,
    -0.0708861117644095;

  Eigen::ArrayXXd rho = get_relabeled_outcomes(observations, num_samples, num_treatments, {});

  double mean_influence_fcn = rho.mean();
  REQUIRE(equal_doubles(mean_influence_fcn, 0, 1e-10));
  REQUIRE(rho.isApprox(rho_expected, 1e-10));
}

TEST_CASE("multi causal relabeling with optional gradient weights works as expected", "[multi causal, relabeling]") {
  size_t num_samples = 10;
  size_t num_treatments = 2;
  std::vector<double> observations = {
    0.747, 0.207, 1.267, -0.503, -1.683, 3.547, 2.237, -5.123,
    -0.493, -0.203, -0.634, -0.304, 1.486, -0.00399999999999999,
    0.056, 1.646, 0.386, -1.344, -0.764, -0.524, 1.012, 0.152, 0.192,
    -0.098, -0.768, 1.582, 0.292, -2.178, 0.492, -0.678
  }; // [Y, W_1, W_2]

  Eigen::ArrayXXd rho = get_relabeled_outcomes(observations, num_samples, num_treatments, {});
  std::vector<double> only_w1 {1, 0};
  Eigen::ArrayXXd rho_only_w1 = get_relabeled_outcomes(observations, num_samples, num_treatments, only_w1);
  std::vector<double> only_w2 {0, 1};
  Eigen::ArrayXXd rho_only_w2 = get_relabeled_outcomes(observations, num_samples, num_treatments, only_w2);
  std::vector<double> times_4 {4, 4};
  Eigen::ArrayXXd rho_times_4 = get_relabeled_outcomes(observations, num_samples, num_treatments, times_4);

  REQUIRE((rho_only_w1.col(0) != 0).all());
  REQUIRE((rho_only_w1.col(1) == 0).all());

  REQUIRE((rho_only_w2.col(0) == 0).all());
  REQUIRE((rho_only_w2.col(1) != 0).all());

  REQUIRE((rho * 4.0 == rho_times_4).all());
}
