#include <map>
#include <unordered_set>
#include <fstream>
#include "utility.h"

#include "catch.hpp"
#include "RelabelingStrategy.h"
#include "InstrumentalRelabelingStrategy.h"
#include "TestUtilities.h"


std::vector<double> get_relabeled_outcomes(Observations observations) {
  std::vector<size_t> sampleIDs;
  for (int i = 0; i < observations.get_num_samples(); i++) {
    sampleIDs.push_back(i);
  }

  std::shared_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy());
  auto relabeled_observations = relabeling_strategy->relabel_outcomes(observations, sampleIDs);

  std::vector<double> relabeled_outcomes;
  for (int i = 0; i < sampleIDs.size(); i++) {
    size_t sampleID = sampleIDs[i];

    REQUIRE(relabeled_observations.count(sampleID));
    relabeled_outcomes.push_back(relabeled_observations.at(sampleID));
  }
  return relabeled_outcomes;
}

TEST_CASE("flipping signs of treatment does not affect relabeled outcomes", "[instrumental, relabeling]") {
  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> flipped_treatment = {0, 1, 1, 1, 0, 1, 0, 1, 1, 1};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};

  Observations observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);
  std::vector<double> first_outcomes = get_relabeled_outcomes(observations);

  Observations flipped_observations = TestUtilities::create_observations(original_outcomes, flipped_treatment, instrument);
  std::vector<double> second_outcomes = get_relabeled_outcomes(flipped_observations);

  REQUIRE(first_outcomes.size() == second_outcomes.size());
  for (int i = 0; i < first_outcomes.size(); i++) {
    double first_outcome = first_outcomes[i];
    double second_outcome = second_outcomes[i];

    REQUIRE(equalDoubles(first_outcome, second_outcome, 1.0e-10));
  }
}

TEST_CASE("scaling instrument scales relabeled outcomes", "[instrumental, relabeling]") {
  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};
  std::vector<double> scaled_instrument = {0, 0, 3, 3, 3, 0, 3, 0, 3, 0};

  Observations observations = TestUtilities::create_observations(original_outcomes,
      treatment, instrument);
  std::vector<double> first_outcomes = get_relabeled_outcomes(observations);

  Observations scaled_observations = TestUtilities::create_observations(original_outcomes,
      treatment, scaled_instrument);
  std::vector<double> second_outcomes = get_relabeled_outcomes(scaled_observations);

  REQUIRE(first_outcomes.size() == second_outcomes.size());
  for (int i = 0; i < first_outcomes.size(); i++) {
    double first_outcome = first_outcomes[i];
    double second_outcome = second_outcomes[i];

    REQUIRE(equalDoubles(3 * first_outcome, second_outcome, 1.0e-10));
  }
}

TEST_CASE("constant treatment leads to no splitting", "[instrumental, relabeling]") {
  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};
  Observations observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);

  std::vector<size_t> sampleIDs;
  for (int i = 0; i < original_outcomes.size(); i++) {
    sampleIDs.push_back(i);
  }

  InstrumentalRelabelingStrategy relabelingStrategy;
  auto relabeled_observations = relabelingStrategy.relabel_outcomes(observations, sampleIDs);

  REQUIRE(relabeled_observations.empty()); // An empty map signals that no splitting should be performed.
}

TEST_CASE("constant instrument leads to no splitting", "[instrumental, relabeling]") {
  std::vector<double> original_outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                           -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  std::vector<double> treatment = {0, 0, 1, 1, 0, 0, 1, 0, 1, 0};
  std::vector<double> instrument = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  Observations observations = TestUtilities::create_observations(original_outcomes, treatment, instrument);

  std::vector<size_t> sampleIDs;
  for (int i = 0; i < original_outcomes.size(); i++) {
    sampleIDs.push_back(i);
  }

  InstrumentalRelabelingStrategy relabelingStrategy;
  auto relabeled_observations = relabelingStrategy.relabel_outcomes(observations, sampleIDs);

  REQUIRE(relabeled_observations.empty()); // An empty map signals that no splitting should be performed.
}
