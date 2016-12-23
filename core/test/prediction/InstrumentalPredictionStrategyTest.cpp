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

  Observations* observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);
  Observations* flipped_observations = TestUtilities::create_observations(original_outcomes,
      flipped_treatment, instrument);

  PredictionStrategy* prediction_strategy = new InstrumentalPredictionStrategy();

  std::vector<double> first_predictions = prediction_strategy->predict(weights_by_sampleID, observations);
  std::vector<double> second_predictions = prediction_strategy->predict(weights_by_sampleID, flipped_observations);

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

  Observations* observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);
  Observations* scaled_observations = TestUtilities::create_observations(original_outcomes,
      treatment, scaled_instrument);

  PredictionStrategy* prediction_strategy = new InstrumentalPredictionStrategy();

  std::vector<double> first_predictions = prediction_strategy->predict(weights_by_sampleID, observations);
  std::vector<double> second_predictions = prediction_strategy->predict(weights_by_sampleID, scaled_observations);

  REQUIRE(first_predictions.size() == 1);
  REQUIRE(second_predictions.size() == 1);
  REQUIRE(equalDoubles(first_predictions[0], second_predictions[0], 1.0e-10));
}