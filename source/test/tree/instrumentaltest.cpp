#include <map>
#include <unordered_set>
#include <fstream>

#include "catch.hpp"
#include "../../src/Tree/RelabelingStrategy.h"
#include "../../src/Tree/InstrumentalRelabelingStrategy.h"

std::vector<double> get_relabeled_outcomes(std::unordered_map<std::string, std::vector<double>>* observations) {
  std::vector<size_t> sampleIDs;
  for (int i = 0; i < (*observations)["outcome"].size(); i++) {
    sampleIDs.push_back(i);
  }

  RelabelingStrategy* relabelingStrategy = new InstrumentalRelabelingStrategy();
  auto relabeled_observations = relabelingStrategy->relabelObservations(observations, sampleIDs);

  std::vector<double> relabeled_outcomes;
  for (int i = 0; i < sampleIDs.size(); i++) {
    size_t sampleID = sampleIDs[i];

    REQUIRE(relabeled_observations.count(sampleID));
    relabeled_outcomes.push_back(relabeled_observations[sampleID]);
  }
  return relabeled_outcomes;
}

TEST_CASE("flipping signs of treatment flips relabeled outcomes", "[instrumental, relabeling]") {
  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> flipped_treatment = {0, 1, 1, 1, 0, 1, 0, 1, 1, 1};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};

  std::unordered_map<std::string, std::vector<double>> observations = {
    {"outcome", original_outcomes}, {"treatment", treatment}, {"instrument", instrument}};
  std::vector<double> outcomes = get_relabeled_outcomes(&observations);

  std::unordered_map<std::string, std::vector<double>> flipped_observations = {
      {"outcome", original_outcomes}, {"treatment", flipped_treatment}, {"instrument", instrument}};
  std::vector<double> flipped_outcomes = get_relabeled_outcomes(&flipped_observations);

  for (int i = 0; i < outcomes.size(); i++) {
    double outcome = outcomes[i];
    double flipped_outcome = flipped_outcomes[i];

    REQUIRE(outcome + flipped_outcome == 0);
  }
}

TEST_CASE("scaling instrument does not affect relabeled outcomes", "[instrumental, relabeling]") {
  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};
  std::vector<double> scaled_instrument = {0, 0, 3, 3, 3, 0, 3, 0, 3, 0};

  std::unordered_map<std::string, std::vector<double>> observations = {
      {"outcome", original_outcomes}, {"treatment", treatment}, {"instrument", instrument}};
  std::vector<double> outcomes = get_relabeled_outcomes(&observations);

  std::unordered_map<std::string, std::vector<double>> scaled_observations = {
      {"outcome", original_outcomes}, {"treatment", treatment}, {"instrument", scaled_instrument}};
  std::vector<double> scaled_outcomes = get_relabeled_outcomes(&scaled_observations);

  for (int i = 0; i < outcomes.size(); i++) {
    double outcome = outcomes[i];
    double scaled_outcome = scaled_outcomes[i];

    REQUIRE(outcome == scaled_outcome);
  }
}
