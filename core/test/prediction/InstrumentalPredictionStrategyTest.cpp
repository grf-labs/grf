/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <map>
#include <unordered_set>
#include <fstream>
#include "Observations.h"
#include "utility.h"
#include "PredictionStrategy.h"
#include "InstrumentalPredictionStrategy.h"
#include "TestUtilities.h"

#include "catch.hpp"

TEST_CASE("flipping signs of treatment flips predictions", "[instrumental, prediction]") {
//  outcomes: {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
//             -5.62082, -9.05911, 3.57729, 3.58593, 8.69386}
//  treatment: {1, 0, 0, 0, 1, 0, 1, 0, 0, 0}
//  instrument: {0, 0, 1, 1, 1, 0, 1, 0, 1, 0}

  std::vector<double> averages = {-1.1251472, 0.3, 0.5, -0.1065444, 0.2};
  std::vector<double> flipped_averages = {-1.1251472, 0.7, 0.5, -0.1065444, 0.3};

  // These values are not required, since we predict using precomputed averages.
  std::unordered_map<size_t, double> weights_by_sampleID;
  Observations observations;

  InstrumentalPredictionStrategy prediction_strategy;
  Prediction first_prediction = prediction_strategy.predict(42, averages, weights_by_sampleID, observations);
  std::vector<double> first_predictions = first_prediction.get_predictions();

  Prediction second_prediction = prediction_strategy.predict(42, flipped_averages, weights_by_sampleID, observations);
  std::vector<double> second_predictions = second_prediction.get_predictions();

  REQUIRE(first_predictions.size() == 1);
  REQUIRE(second_predictions.size() == 1);
  REQUIRE(equal_doubles(first_predictions[0], -second_predictions[0], 1.0e-10));
}
