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

#include "commons/Data.h"
#include "commons/utility.h"
#include "prediction/RegressionPredictionStrategy.h"

#include "catch.hpp"

using namespace grf;

TEST_CASE("flipping signs of outcome flips predictions", "[regression, prediction]") {
  std::vector<double> averages = {1.1251472, 1};
  std::vector<double> flipped_averages = {-1.1251472, 1};

  RegressionPredictionStrategy prediction_strategy;
  std::vector<double> first_prediction = prediction_strategy.predict(averages);
  std::vector<double> second_prediction = prediction_strategy.predict(flipped_averages);

  REQUIRE(first_prediction.size() == 1);
  REQUIRE(second_prediction.size() == 1);
  REQUIRE(equal_doubles(first_prediction[0], -second_prediction[0], 1.0e-10));
}

TEST_CASE("regression variance estimates are positive", "[regression, prediction]") {
  std::vector<double> averages = {1.12, 1};
  std::vector<std::vector<double>> leaf_values = {{3.2, 1}, {4.5, 1}, {6.7, 1}, {-3.5, 1}};

  RegressionPredictionStrategy prediction_strategy;
  std::vector<double> variance = prediction_strategy.compute_variance(
      averages, PredictionValues(leaf_values, 1), 2);

  REQUIRE(variance.size() == 1);
  REQUIRE(variance[0] > 0);
}

TEST_CASE("scaling outcome scales regression variance", "[regression, prediction]") {
  std::vector<double> averages = {2.725, 1};
  std::vector<std::vector<double>> leaf_values = {{3.2, 1}, {4.5, 1}, {6.7, 1}, {-3.5, 1}};

  std::vector<double> scaled_average = {5.45, 1};
  std::vector<std::vector<double>> scaled_leaf_values = {{6.4, 1}, {9.0, 1}, {13.4, 1}, {-7.0, 1}};

  RegressionPredictionStrategy prediction_strategy;
  std::vector<double> first_variance = prediction_strategy.compute_variance(
      averages,
      PredictionValues(leaf_values, 1)
      , 2);
  std::vector<double> second_variance = prediction_strategy.compute_variance(
      scaled_average,
      PredictionValues(scaled_leaf_values, 1), 2);

  REQUIRE(first_variance.size() == 1);
  REQUIRE(second_variance.size() == 1);
  REQUIRE(equal_doubles(first_variance[0], second_variance[0] / 4, 1.0e-10));
}


TEST_CASE("debiased errors are smaller than raw errors", "[regression, prediction]") {
  std::vector<double> average = {2.725, 1};
  std::vector<std::vector<double>> leaf_values = {{3.2, 1}, {4.5, 1}, {6.7, 1}, {-3.5, 1}};

  std::vector<double> outcomes = {6.4, 9.0, 13.4, -7.0};
  Data data(outcomes, 4, 1);
  data.set_outcome_index(0);

  RegressionPredictionStrategy prediction_strategy;

  for (size_t sample=0; sample < 4; ++sample) {
    auto error = prediction_strategy.compute_error(
          sample,
          average,
          PredictionValues(leaf_values, 1),
          data).at(0);
    double debiased_error = error.first;

    // Raw error
    double outcome = data.get_outcome(sample);
    double raw_error = average.at(0) - outcome;
    double mse = raw_error * raw_error;

    REQUIRE(debiased_error < mse);
  }
}
