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
#include "prediction/InstrumentalPredictionStrategy.h"

#include "catch.hpp"

TEST_CASE("flipping signs of treatment flips predictions", "[instrumental, prediction]") {
//  outcomes: {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
//             -5.62082, -9.05911, 3.57729, 3.58593, 8.69386}
//  treatment: {1, 0, 0, 0, 1, 0, 1, 0, 0, 0}
//  instrument: {0, 0, 1, 1, 1, 0, 1, 0, 1, 0}

  std::vector<double> averages = {-1.1251472, 0.3, 0.5, -0.1065444, 0.2};
  std::vector<double> flipped_averages = {-1.1251472, 0.7, 0.5, -0.1065444, 0.3};

  InstrumentalPredictionStrategy prediction_strategy;
  std::vector<double> first_prediction = prediction_strategy.predict(averages);

  std::vector<double> second_prediction = prediction_strategy.predict(flipped_averages);

  REQUIRE(first_prediction.size() == 1);
  REQUIRE(second_prediction.size() == 1);
  REQUIRE(equal_doubles(first_prediction[0], -second_prediction[0], 1.0e-10));
}

TEST_CASE("instrumental variance estimates are positive", "[regression, prediction]") {
  std::vector<double> averages = {1, 0, 4.5, 2, 0.75};
  std::vector<std::vector<double>> leaf_values =
      {{1, 1, 1, 1, 1}, {2, 2, 2, 2, 2}, {-2, -3, 5, -3, -1}, {1, 0, 1, 2, 1}};

  InstrumentalPredictionStrategy prediction_strategy;
  std::vector<double> variance = prediction_strategy.compute_variance(
      averages, PredictionValues(leaf_values, 4, 5), 2);

  REQUIRE(variance.size() == 1);
  REQUIRE(variance[0] > 0);
}

TEST_CASE("scaling outcome scales instrumental variance", "[instrumental, prediction]") {
  std::vector<double> averages = {1, 0, 4.5, 2, 0.75};
  std::vector<std::vector<double>> leaf_values =
      {{1, 1, 1, 1, 1}, {2, 2, 2, 2, 2}, {-2, -3, 5, -3, -1}, {1, 0, 1, 2, 1}};

  std::vector<double> scaled_average = {2, 0, 4.5, 4, 0.75};
  std::vector<std::vector<double>> scaled_leaf_values =
      {{2, 1, 1, 2, 1}, {4, 2, 2, 4, 2}, {-4, -3, 5, -6, -1}, {2, 0, 1, 4, 1}};

  InstrumentalPredictionStrategy prediction_strategy;
  std::vector<double> first_variance = prediction_strategy.compute_variance(
      averages,
      PredictionValues(leaf_values, 4, 5),
      2);
  std::vector<double> second_variance = prediction_strategy.compute_variance(
      scaled_average,
      PredictionValues(scaled_leaf_values, 4, 5),
      2);

  REQUIRE(first_variance.size() == 1);
  REQUIRE(second_variance.size() == 1);
  REQUIRE(equal_doubles(first_variance[0], second_variance[0] / 4, 10e-10));
}
