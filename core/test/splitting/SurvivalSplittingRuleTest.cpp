/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include "splitting/SurvivalSplittingRule.h"
#include "relabeling/NoopRelabelingStrategy.h"

#include "commons/utility.h"
#include "utilities/ForestTestUtilities.h"
#include "utilities/FileTestUtilities.h"

#include "catch.hpp"

using namespace grf;

// Splitting rule input setup to emulate one run of node zero (all data) splitting on all features
// returning a vector of the best logrank statistic for each feature.
std::vector<double> run_splits(const Data& data,
                               const TreeOptions& options,
                               const std::unique_ptr<SurvivalSplittingRule>& splitting_rule,
                               const std::unique_ptr<RelabelingStrategy>& relabeling_strategy,
                               size_t num_features) {
  std::vector<double> best_logranks;

  size_t node = 0;
  size_t size_node = data.get_num_rows();
  Eigen::ArrayXXd responses_by_sample(size_node, 1);
  std::vector<std::vector<size_t>> samples(1);
  for (size_t sample = 0; sample < size_node; ++sample) {
    samples[node].push_back(sample);
  }
  relabeling_strategy->relabel(samples[node], data, responses_by_sample);
  double split_value = 0;
  size_t split_variable = 0;
  bool send_missing_left = true;

  for (size_t split_var = 0; split_var < num_features; split_var++) {
    std::vector<size_t> possible_split_vars;
    possible_split_vars.push_back(split_var);
    double best_logrank = 0;
    splitting_rule->find_best_split_internal(data,
                                             possible_split_vars,
                                             responses_by_sample,
                                             samples[node],
                                             split_value,
                                             split_variable,
                                             send_missing_left,
                                             best_logrank);
    best_logranks.push_back(best_logrank);
   }

   return best_logranks;
}

TEST_CASE("survival splitting logrank calculation is correct", "[survival], [splitting]") {
  auto data_vec = load_data("test/splitting/resources/survival_data_logrank.csv");
  Data data(data_vec);
  size_t num_features = 500;
  data.set_outcome_index(num_features);
  data.set_censor_index(num_features + 1);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  std::unique_ptr<SurvivalSplittingRule> surv_splitting_rule(new SurvivalSplittingRule(options.get_alpha()));

  std::vector<double> logranks = run_splits(data, options, surv_splitting_rule, relabeling_strategy, num_features);

  std::vector<std::vector<double>> expected_logranks = FileTestUtilities::read_csv_file(
      "test/splitting/resources/survival_data_logrank_expected.csv");

  for (size_t i = 0; i < logranks.size(); ++i) {
    double logrank = logranks[i];
    double expected_logrank = expected_logranks[i][0];
    REQUIRE(equal_doubles(logrank, expected_logrank, 1e-6));
  }
}
