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
#include "relabeling/RelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"

using namespace grf;

TEST_CASE("simple quantile relabeling", "[quantile, relabeling]") {
  std::vector<double> outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                  -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  Data data(outcomes, 10, 1);
  data.set_outcome_index(0);

  std::vector<size_t> samples;
  for (size_t i = 0; i < data.get_num_rows(); ++i) {
    samples.push_back(i);
  }

  QuantileRelabelingStrategy relabeling_strategy({0.25, 0.5, 0.75});

  Eigen::ArrayXXd relabeled_observations(data.get_num_rows(), 1);
  bool stop = relabeling_strategy.relabel(samples, data, relabeled_observations);
  REQUIRE(stop == false);

  std::vector<double> relabeled_outcomes;
  for (auto& sample : samples) {
    relabeled_outcomes.push_back(relabeled_observations(sample));
  }

  std::vector<double> expected_outcomes = {0, 0, 3, 1, 2, 1, 0, 2, 2, 3};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}

TEST_CASE("quantile relabeling subset of observations", "[quantile, relabeling]") {
  std::vector<double> outcomes = {-2.32996, 0.388327, 6.61931, -9.30856, -8.93077,
                                  0.594004, 3.42299, -9.84604, -2.33169, -8.66316};
  Data data(outcomes, 10, 1);
  data.set_outcome_index(0);

  std::vector<size_t> samples = {1, 3, 5, 7, 9};

  QuantileRelabelingStrategy relabeling_strategy({0.5, 0.75});

  Eigen::ArrayXXd relabeled_observations(data.get_num_rows(), 1);
  bool stop = relabeling_strategy.relabel(samples, data, relabeled_observations);
  REQUIRE(stop == false);

  std::vector<double> relabeled_outcomes;
  for (auto& sample : samples) {
    relabeled_outcomes.push_back(relabeled_observations(sample));
  }

  std::vector<double> expected_outcomes = {1, 0, 2, 0, 0};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}
