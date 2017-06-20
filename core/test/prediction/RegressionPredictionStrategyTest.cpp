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

#include <map>
#include <unordered_set>
#include <fstream>
#include "commons/Observations.h"
#include "commons/utility.h"
#include "prediction/RegressionPredictionStrategy.h"

#include "catch.hpp"

TEST_CASE("flipping signs of outcome flips predictions", "[regression, prediction]") {
  std::vector<double> averages = {1.1251472};
  std::vector<double> flipped_averages = {-1.1251472};

  RegressionPredictionStrategy prediction_strategy;
  std::vector<double> first_prediction = prediction_strategy.predict(averages);
  std::vector<double> second_prediction = prediction_strategy.predict(flipped_averages);

  REQUIRE(first_prediction.size() == 1);
  REQUIRE(second_prediction.size() == 1);
  REQUIRE(equal_doubles(first_prediction[0], -second_prediction[0], 1.0e-10));
}

TEST_CASE("regression variance estimates are positive", "[regression, prediction]") {
  std::vector<double> averages = {1.12};
  std::vector<std::vector<double>> leaf_values = {{3.2}, {4.5}, {6.7}, {-3.5}};

  RegressionPredictionStrategy prediction_strategy;
  std::vector<double> variance = prediction_strategy.compute_variance(
      averages, PredictionValues(leaf_values, 4, 1), 2);

  REQUIRE(variance.size() == 1);
  REQUIRE(variance[0] > 0);
}

TEST_CASE("scaling outcome scales regression variance", "[regression, prediction]") {
  std::vector<double> averages = {2.725};
  std::vector<std::vector<double>> leaf_values = {{3.2}, {4.5}, {6.7}, {-3.5}};

  std::vector<double> scaled_average = {5.45};
  std::vector<std::vector<double>> scaled_leaf_values = {{6.4}, {9.0}, {13.4}, {-7.0}};

  RegressionPredictionStrategy prediction_strategy;
  std::vector<double> first_variance = prediction_strategy.compute_variance(
      averages,
      PredictionValues(leaf_values, 4, 1)
      , 2);
  std::vector<double> second_variance = prediction_strategy.compute_variance(
      scaled_average,
      PredictionValues(scaled_leaf_values, 4, 1), 2);

  REQUIRE(first_variance.size() == 1);
  REQUIRE(second_variance.size() == 1);
  REQUIRE(equal_doubles(first_variance[0], second_variance[0] / 4, 1.0e-10));
}

