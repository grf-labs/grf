/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <map>
#include <unordered_set>
#include <fstream>

#include "catch.hpp"
#include "relabeling/RelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "utilities/TestUtilities.h"

TEST_CASE("simple quantile relabeling", "[quantile, relabeling]") {
  std::vector<double> outcomes = {-9.99984, -7.36924, 5.11211, -0.826997, 0.655345,
                                  -5.62082, -9.05911, 3.57729, 3.58593, 8.69386};
  Observations observations = TestUtilities::create_observations(outcomes);

  std::vector<size_t> sampleIDs;
  for (size_t i = 0; i < outcomes.size(); ++i) {
    sampleIDs.push_back(i);
  }

  QuantileRelabelingStrategy relabeling_strategy({0.25, 0.5, 0.75});
  auto relabeled_observations = relabeling_strategy.relabel(sampleIDs, observations);

  std::vector<double> relabeled_outcomes;
  for (auto& sampleID : sampleIDs) {
    REQUIRE(relabeled_observations.count(sampleID));
    relabeled_outcomes.push_back(relabeled_observations.at(sampleID));
  }

  std::vector<double> expected_outcomes = {0, 0, 3, 1, 2, 1, 0, 2, 2, 3};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}

TEST_CASE("quantile relabeling subset of observations", "[quantile, relabeling]") {
  std::vector<double> outcomes = {-2.32996, 0.388327, 6.61931, -9.30856, -8.93077,
                                  0.594004, 3.42299, -9.84604, -2.33169, -8.66316};
  Observations observations = TestUtilities::create_observations(outcomes);

  std::vector<size_t> sampleIDs = {1, 3, 5, 7, 9};

  QuantileRelabelingStrategy relabeling_strategy({0.5, 0.75});
  auto relabeled_observations = relabeling_strategy.relabel(sampleIDs, observations);

  std::vector<double> relabeled_outcomes;
  for (auto& sampleID : sampleIDs) {
    REQUIRE(relabeled_observations.count(sampleID));
    relabeled_outcomes.push_back(relabeled_observations.at(sampleID));
  }

  std::vector<double> expected_outcomes = {1, 0, 2, 0, 0};
  REQUIRE(expected_outcomes == relabeled_outcomes);
}
