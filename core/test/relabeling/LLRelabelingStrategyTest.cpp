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

#include <map>
#include <unordered_set>
#include <fstream>

#include "commons/utility.h"
#include "commons/DefaultData.h"
#include "catch.hpp"
#include "relabeling/RelabelingStrategy.h"
#include "relabeling/LLRelabelingStrategy.h"

using namespace grf;

TEST_CASE("simple local linear relabeling", "[ll_regression, relabeling]") {
  std::unique_ptr<Data> data = load_data("test/forest/resources/friedman.csv");
  std::vector<std::vector<double>> expected_outcomes = FileTestUtilities::read_csv_file(
          "test/forest/resources/ll_relabeled_outcomes.csv");

  data->set_outcome_index(10);
  std::vector<size_t> ll_split_variables = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> split_lambda = {0.1};
  bool weight_penalty = false;

  std::vector<size_t> samples;
  for (size_t i = 0; i < data.get_num_rows(); ++i) {
    samples.push_back(i);
  }

  LLRelabelingStrategy relabeling_strategy(split_lambda, weight_penalty, ll_split_variables);

  std::vector<double> relabeled_observations(data.get_num_rows());
  bool stop = relabeling_strategy.relabel(samples, data, relabeled_observations);
  REQUIRE(stop == false);

  std::vector<double> relabeled_outcomes;
  for (auto& sample : samples) {
    relabeled_outcomes.push_back(relabeled_observations.at(sample));
  }

  REQUIRE(expected_outcomes == relabeled_outcomes);
}
