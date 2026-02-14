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

#include "prediction/DefaultPredictionStrategy.h"
#include "prediction/QuantilePredictionStrategy.h"


#include "catch.hpp"

using namespace grf;

TEST_CASE("simple quantile prediction", "[quantile, prediction]") {
  std::pair<std::vector<size_t>, std::vector<double>> weights_by_sample = {
      {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
       {0.0, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.0, 0.1}}};

  std::vector<double> outcomes = { -9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                   -5.62082, -9.05911, 3.57729, 3.58593, 8.69386 };
  Data data(outcomes, 10, 1);
  data.set_outcome_index(0);

  QuantilePredictionStrategy prediction_strategy({0.25, 0.5, 0.75});
  std::vector<double> predictions =  prediction_strategy.predict(0, weights_by_sample, data, data);

  std::vector<double> expected_predictions = {-7.36924, -0.826997, 5.11211};
  REQUIRE(predictions == expected_predictions);
}

TEST_CASE("prediction with skewed quantiles", "[quantile, prediction]") {
  std::pair<std::vector<size_t>, std::vector<double>> weights_by_sample = {
      {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
       {0.0, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.0, 0.1}}};

  std::vector<double> outcomes =  { -1.99984, -0.36924, 0.11211, -1.826997, 1.655345,
                                    -1.62082, -0.05911, 0.57729, 0.58593, 1.69386 };
  Data data(outcomes, 10, 1);
  data.set_outcome_index(0);

  QuantilePredictionStrategy prediction_strategy({0.5, 0.75, 0.80, 0.90});
  std::vector<double> predictions = prediction_strategy.predict(42, weights_by_sample, data, data);

  // Check that all predictions fall within a reasonable range.
  for (auto& prediction : predictions) {
    REQUIRE(-2.0 < prediction);
    REQUIRE(prediction < 2.0);
  }
}

TEST_CASE("prediction with repeated quantiles", "[quantile, prediction]") {
  std::pair<std::vector<size_t>, std::vector<double>> weights_by_sample = {
      {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
       {0.0, 0.1, 0.2, 0.1, 0.1, 0.1, 0.2, 0.1, 0.0, 0.1}}};

  std::vector<double> outcomes  = { -9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                    -5.62082, -9.05911, 3.57729, 3.58593, 8.69386 };
  Data data(outcomes, 10, 1);
  data.set_outcome_index(0);

  std::vector<double> first_predictions = QuantilePredictionStrategy({0.5})
          .predict(42, weights_by_sample, data, data);
  std::vector<double> second_predictions = QuantilePredictionStrategy({0.25, 0.5, 0.75})
      .predict(42, weights_by_sample, data, data);
  std::vector<double> third_predictions = QuantilePredictionStrategy({0.5, 0.5, 0.5})
      .predict(42, weights_by_sample, data, data);

  REQUIRE(first_predictions[0] == second_predictions[1]);
  for (auto prediction : third_predictions) {
    REQUIRE(prediction == first_predictions[0]);
  }
}
