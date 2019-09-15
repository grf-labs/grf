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

#include <fstream>
#include <map>
#include <memory>
#include <unordered_set>

#include "catch.hpp"
#include "commons/utility.h"
#include "relabeling/RelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"

using namespace grf;

std::vector<double> get_relabeled_outcomes(std::vector<double> observations, size_t num_samples) {
  DefaultData data(observations, num_samples, 3);
  data.set_outcome_index(0);
  data.set_treatment_index(1);
  data.set_instrument_index(2);

  std::vector<size_t> samples;
  for (size_t i = 0; i < num_samples; ++i) {
    samples.push_back(i);
  }

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy());

  std::vector<double> relabeled_observations(num_samples);
  bool stop = relabeling_strategy->relabel(samples, data, relabeled_observations);
  if (stop) {
    return std::vector<double>();
  }

  std::vector<double> relabeled_outcomes;
  relabeled_outcomes.reserve(samples.size());
  for (auto& sample : samples) {
    relabeled_outcomes.push_back(relabeled_observations.at(sample));
  }
  return relabeled_outcomes;
}

TEST_CASE("flipping signs of treatment does not affect relabeled outcomes", "[instrumental, relabeling]") {
  std::vector<double> observations = {
      -9.99984, -7.36924, 5.11211, -0.826997, 0.655345, -5.62082, -9.05911, 3.57729, 3.58593, 8.69386, // outcomes
      1, 0, 0, 0, 1, 0, 1, 0, 0, 0, // treatment
      0, 0, 1, 1, 1, 0, 1, 0, 1, 0 }; // instrument

  std::vector<double> flipped_observations = {
      -9.99984, -7.36924, 5.11211, -0.826997, 0.655345, -5.62082, -9.05911, 3.57729, 3.58593, 8.69386, // outcomes
      0, 1, 1, 1, 0, 1, 0, 1, 1, 1, // treatment
      0, 0, 1, 1, 1, 0, 1, 0, 1, 0 }; // instrument

  std::vector<double> first_outcomes = get_relabeled_outcomes(observations, 10);
  std::vector<double> second_outcomes = get_relabeled_outcomes(flipped_observations, 10);

  REQUIRE(first_outcomes.size() == second_outcomes.size());
  for (size_t i = 0; i < first_outcomes.size(); ++i) {
    double first_outcome = first_outcomes[i];
    double second_outcome = second_outcomes[i];

    REQUIRE(equal_doubles(first_outcome, second_outcome, 1.0e-10));
  }
}

TEST_CASE("scaling instrument scales relabeled outcomes", "[instrumental, relabeling]") {
  std::vector<double> outcomes = { };
  std::vector<double> treatment = {1, 0, 0, 0, 1, 0, 1, 0, 0, 0};
  std::vector<double> instrument = {0, 0, 1, 1, 1, 0, 1, 0, 1, 0};
  std::vector<double> scaled_instrument = {0, 0, 3, 3, 3, 0, 3, 0, 3, 0};

  std::vector<double> observations = {
      -9.99984, -7.36924, 5.11211, -0.826997, 0.655345, -5.62082, -9.05911, 3.57729, 3.58593, 8.69386, // outcomes
      1, 0, 0, 0, 1, 0, 1, 0, 0, 0, // treatment
      0, 0, 1, 1, 1, 0, 1, 0, 1, 0 }; // instrument

  std::vector<double> scaled_observations = {
      -9.99984, -7.36924, 5.11211, -0.826997, 0.655345, -5.62082, -9.05911, 3.57729, 3.58593, 8.69386, // outcomes
      1, 0, 0, 0, 1, 0, 1, 0, 0, 0, // treatment
      0, 0, 3, 3, 3, 0, 3, 0, 3, 0 }; // scaled instrument

  std::vector<double> first_outcomes = get_relabeled_outcomes(observations, 10);
  std::vector<double> second_outcomes = get_relabeled_outcomes(scaled_observations, 10);

  REQUIRE(first_outcomes.size() == second_outcomes.size());
  for (size_t i = 0; i < first_outcomes.size(); ++i) {
    double first_outcome = first_outcomes[i];
    double second_outcome = second_outcomes[i];

    REQUIRE(equal_doubles(3 * first_outcome, second_outcome, 1.0e-10));
  }
}

TEST_CASE("constant treatment leads to no splitting", "[instrumental, relabeling]") {
  std::vector<double> observations = {
      -9.99984, -7.36924, 5.11211, -0.826997, 0.655345, -5.62082, -9.05911, 3.57729, 3.58593, 8.69386, // outcomes
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // treatment
      0, 0, 1, 1, 1, 0, 1, 0, 1, 0 }; // instrument

  std::vector<double> relabeled_outcomes = get_relabeled_outcomes(observations, 10);
  REQUIRE(relabeled_outcomes.empty()); // An empty map signals that no splitting should be performed.
}

TEST_CASE("constant instrument leads to no splitting", "[instrumental, relabeling]") {
  std::vector<double> observations = {
      -9.99984, -7.36924, 5.11211, -0.826997, 0.655345, -5.62082, -9.05911, 3.57729, 3.58593, 8.69386, // outcomes
      0, 0, 1, 1, 0, 0, 1, 0, 1, 0, // treatment
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }; // instrument

  std::vector<double> relabeled_outcomes = get_relabeled_outcomes(observations, 10);
  REQUIRE(relabeled_outcomes.empty()); // An empty map signals that no splitting should be performed.
}
