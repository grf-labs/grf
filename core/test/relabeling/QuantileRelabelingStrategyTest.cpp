#include <map>
#include <unordered_set>
#include <fstream>

#include "catch.hpp"
#include "RelabelingStrategy.h"
#include "QuantileRelabelingStrategy.h"

TEST_CASE("simple quantile relabeling", "[quantile, relabeling]") {
  std::vector<double> outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                  -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::unordered_map<std::string, std::vector<double>> observations = {{"outcome", outcomes}};

  std::vector<size_t> sampleIDs;
  for (int i = 0; i < outcomes.size(); i++) {
    sampleIDs.push_back(i);
  }

  std::vector<double>* quantiles = new std::vector<double>({0.25, 0.5, 0.75});
  RelabelingStrategy* relabelingStrategy = new QuantileRelabelingStrategy(quantiles);
  auto relabeled_observations = relabelingStrategy->relabelObservations(&observations, sampleIDs);

  std::vector<double> relabeled_outcomes;
  for (int i = 0; i < sampleIDs.size(); i++) {
    size_t sampleID = sampleIDs[i];

    REQUIRE(relabeled_observations.count(sampleID));
    relabeled_outcomes.push_back(relabeled_observations[sampleID]);
  }

  std::vector<double> expected_outcomes = {0, 0, 3, 1, 2, 1, 0, 2, 2, 3};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}

TEST_CASE("quantile relabeling subset of observations", "[quantile, relabeling]") {
  std::vector<double> outcomes = {-2.32996, 0.388327, 6.61931, -9.30856, -8.93077,
                                  0.594004, 3.42299, -9.84604, -2.33169, -8.66316};
  std::unordered_map<std::string, std::vector<double>> observations = {{"outcome", outcomes}};

  std::vector<size_t> sampleIDs = {1, 3, 5, 7, 9};

  std::vector<double>* quantiles = new std::vector<double>({0.5, 0.75});
  RelabelingStrategy* relabelingStrategy = new QuantileRelabelingStrategy(quantiles);
  auto relabeled_observations = relabelingStrategy->relabelObservations(&observations, sampleIDs);

  std::vector<double> relabeled_outcomes;
  for (int i = 0; i < sampleIDs.size(); i++) {
    size_t sampleID = sampleIDs[i];

    REQUIRE(relabeled_observations.count(sampleID));
    relabeled_outcomes.push_back(relabeled_observations[sampleID]);
  }

  std::vector<double> expected_outcomes = {1, 0, 2, 0, 0};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}
