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

#include "splitting/factory/SurvivalSplittingRuleFactory.h"
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
                               const std::unique_ptr<SplittingRuleFactory>& splitting_rule_factory,
                               const std::unique_ptr<RelabelingStrategy>& relabeling_strategy,
                               size_t num_features) {
  std::vector<double> best_logranks;

  std::unique_ptr<SplittingRule> splitting_rule = splitting_rule_factory->create(data.get_num_rows(), options);
  size_t node = 0;
  size_t size_node = data.get_num_rows();
  std::vector<double> responses_by_sample(size_node);
  std::vector<std::vector<size_t>> samples(1);
  for (size_t sample = 0; sample < size_node; ++sample) {
    samples[node].push_back(sample);
  }
  relabeling_strategy->relabel(samples[node], data, responses_by_sample);
  std::vector<size_t> split_vars(1);
  std::vector<double> split_values(1);
  std::vector<bool> send_missing_left(1);

  for (size_t split_var = 0; split_var < num_features; split_var++) {
    std::vector<size_t> possible_split_vars;
    possible_split_vars.push_back(split_var);

    splitting_rule->find_best_split(data,
                                    node,
                                    possible_split_vars,
                                    responses_by_sample,
                                    samples,
                                    split_vars,
                                    split_values,
                                    send_missing_left);

    best_logranks.push_back(splitting_rule->test_statistic);
   }

   return best_logranks;
}

TEST_CASE("survival splitting logrank calculation is correct", "[survival], [splitting]") {
  std::unique_ptr<Data> data = load_data("test/splitting/resources/survival_data_logrank.csv");
  size_t num_features = 500;
  data->set_outcome_index(num_features);
  data->set_censor_index(num_features + 1);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new SurvivalSplittingRuleFactory());

  std::vector<double> logranks = run_splits(*data, options, splitting_rule_factory, relabeling_strategy, num_features);

  std::vector<std::vector<double>> expected_logranks = FileTestUtilities::read_csv_file(
      "test/splitting/resources/survival_data_logrank_expected.csv");

  for (size_t i = 0; i < logranks.size(); ++i) {
    double logrank = logranks[i];
    double expected_logrank = expected_logranks[i][0];
    REQUIRE(equal_doubles(logrank, expected_logrank, 1e-6));
  }
}
