#include <map>
#include <unordered_set>
#include <fstream>

#include "catch.hpp"
#include "../../src/Tree/RelabelingStrategy.h"
#include "../../src/Tree/QuantileRelabelingStrategy.h"

double fRand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void printRand() {
  for (int i = 0; i < 10; i++) {
    std::cout << fRand(-10.0, 10.0) << ", ";
  }
  std::cout << std::endl;
}

TEST_CASE("simple quantile relabeling succeeds", "[quantile, relabeling]") {
  std::vector<double>* quantiles = new std::vector<double>({0.25, 0.5, 0.75});
  RelabelingStrategy* relabelingStrategy = new QuantileRelabelingStrategy(quantiles);

std::vector<double> outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::unordered_map<std::string, std::vector<double>> observations = {{"outcome", outcomes}};

  std::vector<size_t> sampleIDs;
  for (int i = 0; i < outcomes.size(); i++) {
    sampleIDs.push_back(i);
  }

  auto relabeled_observations = relabelingStrategy->relabelObservations(&observations, sampleIDs);

  std::vector<double> relabeled_outcomes;
  for (int i = 0; i < outcomes.size(); i++) {
    REQUIRE(relabeled_observations.count(i));
    relabeled_outcomes.push_back(relabeled_observations[i]);
  }

  std::vector<double> expected_outcomes = {0, 0, 3, 1, 2, 1, 0, 2, 2, 3};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}

