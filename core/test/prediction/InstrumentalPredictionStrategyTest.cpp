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

  Author: Julie Tibshirani (jtibs@cs.stanford.edu)
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
  std::unordered_map<size_t, double> weights_by_sampleID = {
      {0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.1}, {4, 0.1},
      {5, 0.1}, {6, 0.2}, {7, 0.1}, {8, 0.0}, {9, 0.1}};

  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> flipped_treatment = {0, 1, 1, 1, 0, 1, 0, 1, 1, 1};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};

  Observations observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);
  Observations flipped_observations = TestUtilities::create_observations(original_outcomes,
      flipped_treatment, instrument);

  InstrumentalPredictionStrategy prediction_strategy;
  Prediction first_prediction = prediction_strategy.predict({}, weights_by_sampleID, observations);
  std::vector<double> first_predictions = first_prediction.get_predictions();

  Prediction second_prediction = prediction_strategy.predict({}, weights_by_sampleID, flipped_observations);
  std::vector<double> second_predictions = second_prediction.get_predictions();

  REQUIRE(first_predictions.size() == 1);
  REQUIRE(second_predictions.size() == 1);
  REQUIRE(equalDoubles(first_predictions[0], -second_predictions[0], 1.0e-10));
}

TEST_CASE("scaling instrument does not affect prediction", "[instrumental, prediction]") {
  std::unordered_map<size_t, double> weights_by_sampleID = {
      {0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.1}, {4, 0.1},
      {5, 0.1}, {6, 0.2}, {7, 0.1}, {8, 0.0}, {9, 0.1}};

  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};
  std::vector<double> scaled_instrument = {0, 0, 3, 3, 3, 0, 3, 0, 3, 0};

  Observations observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);
  Observations scaled_observations = TestUtilities::create_observations(original_outcomes,
      treatment, scaled_instrument);

  InstrumentalPredictionStrategy prediction_strategy;
  Prediction first_prediction = prediction_strategy.predict({}, weights_by_sampleID, observations);
  std::vector<double> first_predictions = first_prediction.get_predictions();

  Prediction second_prediction = prediction_strategy.predict({}, weights_by_sampleID, scaled_observations);
  std::vector<double> second_predictions = second_prediction.get_predictions();

  REQUIRE(first_predictions.size() == 1);
  REQUIRE(second_predictions.size() == 1);
  REQUIRE(equalDoubles(first_predictions[0], second_predictions[0], 1.0e-10));
}

TEST_CASE("predicting with variance gives same predictions", "[instrumental, prediction]") {
  std::vector<std::vector<size_t>> leaf_sampleIDs = {
      {1, 2, 4, 5}, {}, {2, 5}, {0}, {0, 1}, {5, 6, 7, 8}};

  std::unordered_map<size_t, double> weights_by_sampleID = {
      {0, 6/20.0}, {1, 3/20.0}, {2, 3/20.0}, {3, 0.0}, {4, 1/20.0},
      {5, 4/20.0}, {6, 1/20.0}, {7, 1/20.0}, {8, 1/20.0}, {9, 0.0}};

  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};

  Observations observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);

  InstrumentalPredictionStrategy prediction_strategy;
  Prediction first_prediction = prediction_strategy.predict({}, weights_by_sampleID, observations);
  std::vector<double> first_predictions = first_prediction.get_predictions();

  Prediction second_prediction = prediction_strategy.predict_with_variance(leaf_sampleIDs, observations, 2);
  std::vector<double> second_predictions = second_prediction.get_predictions();

  REQUIRE(first_predictions.size() == 1);
  REQUIRE(second_predictions.size() == 1);
  REQUIRE(equalDoubles(first_predictions[0], second_predictions[0], 1.0e-10));
}