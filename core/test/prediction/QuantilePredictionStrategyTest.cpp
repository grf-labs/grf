#include <map>
#include <unordered_set>
#include <fstream>
#include "PredictionStrategy.h"
#include "QuantilePredictionStrategy.h"
#include "TestUtilities.h"

#include "catch.hpp"

TEST_CASE("simple quantile prediction", "[quantile, prediction]") {
  std::unordered_map<size_t, double> weights_by_sampleID = {
      {0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.1}, {4, 0.1},
      {5, 0.1}, {6, 0.2}, {7, 0.1}, {8, 0.0}, {9, 0.1}};

  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  Observations observations = TestUtilities::create_observations(original_outcomes);

  QuantilePredictionStrategy prediction_strategy({0.25, 0.5, 0.75});
  std::vector<double> predictions = prediction_strategy.predict({}, weights_by_sampleID, observations);

  std::vector<double> expected_predictions = {-7.36924, -0.826997, 5.11211};
  REQUIRE(predictions == expected_predictions);
}

TEST_CASE("prediction with skewed quantiles", "[quantile, prediction]") {
  std::unordered_map<size_t, double> weights_by_sampleID = {
      {0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.1}, {4, 0.1},
      {5, 0.1}, {6, 0.2}, {7, 0.1}, {8, 0.0}, {9, 0.1}};

  std::vector<double> original_outcomes = {-1.99984, -0.36924, 0.11211, -1.826997, 1.655345,
                                           -1.62082, -0.05911, 0.57729, 0.58593, 1.69386};
  Observations observations = TestUtilities::create_observations(original_outcomes);

  QuantilePredictionStrategy prediction_strategy({0.5, 0.75, 0.80, 0.90});
  std::vector<double> predictions = prediction_strategy.predict({}, weights_by_sampleID, observations);

  // Check that all predictions fall within a reasonable range.
  for (auto &prediction : predictions) {
    REQUIRE(-2.0 < prediction);
    REQUIRE(prediction < 2.0);
  }
}

TEST_CASE("prediction with repeated quantiles", "[quantile, prediction]") {
  std::unordered_map<size_t, double> weights_by_sampleID = {
      {0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.1}, {4, 0.1},
      {5, 0.1}, {6, 0.2}, {7, 0.1}, {8, 0.0}, {9, 0.1}};

  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  Observations observations = TestUtilities::create_observations(original_outcomes);

  std::vector<double> first_predictions = QuantilePredictionStrategy({0.5})
          .predict({}, weights_by_sampleID, observations);
  std::vector<double> second_predictions = QuantilePredictionStrategy({0.25, 0.5, 0.75})
      .predict({}, weights_by_sampleID, observations);
  std::vector<double> third_predictions = QuantilePredictionStrategy({0.5, 0.5, 0.5})
      .predict({}, weights_by_sampleID, observations);

  REQUIRE(first_predictions[0] == second_predictions[1]);
  for (auto prediction : third_predictions) {
    REQUIRE(prediction == first_predictions[0]);
  }
}