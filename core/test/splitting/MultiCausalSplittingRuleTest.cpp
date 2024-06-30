/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#include "splitting/SplittingRule.h"
#include "splitting/InstrumentalSplittingRule.h"
#include "splitting/MultiCausalSplittingRule.h"
#include "relabeling/NoopRelabelingStrategy.h"

#include "commons/utility.h"
#include "utilities/ForestTestUtilities.h"
#include "utilities/FileTestUtilities.h"

#include "catch.hpp"

using namespace grf;

// Splitting rule input setup to emulate one run of node zero (first `size_node` observations) splitting
// on all features returning a vector containing the best split variable, best split value, and missing direction.
std::vector<double> run_splits_multi(const Data& data,
                                     size_t size_node,
                                     const TreeOptions& options,
                                     const std::unique_ptr<SplittingRule>& splitting_rule,
                                     const std::unique_ptr<RelabelingStrategy>& relabeling_strategy,
                                     size_t num_features) {
  size_t node = 0;
  Eigen::ArrayXXd responses_by_sample(size_node, data.get_num_outcomes());
  std::vector<std::vector<size_t>> samples(1);
  for (size_t sample = 0; sample < size_node; ++sample) {
    samples[node].push_back(sample);
  }
  relabeling_strategy->relabel(samples[node], data, responses_by_sample);

  std::vector<size_t> possible_split_vars;
  for (size_t j = 0; j < num_features; j++) {
    possible_split_vars.push_back(j);
  }
  std::vector<size_t> split_vars(1);
  std::vector<double> split_values(1);
  std::vector<bool> send_missing_left(1);

  splitting_rule->find_best_split(data,
                                  node,
                                  possible_split_vars,
                                  responses_by_sample,
                                  samples,
                                  split_vars,
                                  split_values,
                                  send_missing_left);

return {(double) split_vars[0], split_values[0], (double) send_missing_left[0]};
}

TEST_CASE("multi causal splitting with one treatment is identical to instrumental splitting", "[multi_causal], [splitting]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);
  size_t num_features = 10;
  size_t num_treatments = data.get_num_treatments();

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  // Ignore relabeling for the purpose on these tests and only split directly on outcomes.
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto inst_splitting_rule = std::unique_ptr<SplittingRule>(new InstrumentalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty()));
  auto multi_splitting_rule = std::unique_ptr<SplittingRule>(new MultiCausalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      1,
      num_treatments));


  for (auto size_node : {1000, 500, 250, 113, 25, 10, 3}) {
    std::vector<double> inst = run_splits_multi(data, size_node, options, inst_splitting_rule, relabeling_strategy, num_features);
    std::vector<double> multi = run_splits_multi(data, size_node, options, multi_splitting_rule, relabeling_strategy, num_features);

    REQUIRE(equal_doubles(inst[0], multi[0], 1e-6));
    REQUIRE(equal_doubles(inst[1], multi[1], 1e-6));
    REQUIRE(equal_doubles(inst[2], multi[2], 1e-6));
  }
}

TEST_CASE("multi causal splitting with one treatment on NaN data is identical to instrumental splitting", "[multi_causal], [splitting], [NaN]") {
  auto data_vec = load_data("test/forest/resources/causal_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);
  size_t num_features = 10;
  size_t num_treatments = data.get_num_treatments();

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto inst_splitting_rule = std::unique_ptr<SplittingRule>(new InstrumentalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty()));
  auto multi_splitting_rule = std::unique_ptr<SplittingRule>(new MultiCausalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      1,
      num_treatments));

  for (auto size_node : {1000, 500, 250, 113, 25, 10, 3}) {
    std::vector<double> inst = run_splits_multi(data, size_node, options, inst_splitting_rule, relabeling_strategy, num_features);
    std::vector<double> multi = run_splits_multi(data, size_node, options, multi_splitting_rule, relabeling_strategy, num_features);

    REQUIRE(equal_doubles(inst[0], multi[0], 1e-6));
    REQUIRE(equal_doubles(inst[1], multi[1], 1e-6));
    REQUIRE(equal_doubles(inst[2], multi[2], 1e-6));
  }
}

TEST_CASE("multi causal splitting with many identical treatments is identical to instrumental splitting", "[multi_causal], [splitting]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index({11, 11, 11});
  data.set_instrument_index(11);
  size_t num_features = 10;
  size_t num_treatments = data.get_num_treatments();

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  // Ignore relabeling for the purpose on these tests and only split directly on outcomes.
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto inst_splitting_rule = std::unique_ptr<SplittingRule>(new InstrumentalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty()));
  auto multi_splitting_rule = std::unique_ptr<SplittingRule>(new MultiCausalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      1,
      num_treatments));


  for (auto size_node : {1000, 500, 250, 113, 25, 10, 3}) {
    std::vector<double> inst = run_splits_multi(data, size_node, options, inst_splitting_rule, relabeling_strategy, num_features);
    std::vector<double> multi = run_splits_multi(data, size_node, options, multi_splitting_rule, relabeling_strategy, num_features);

    REQUIRE(equal_doubles(inst[0], multi[0], 1e-6));
    REQUIRE(equal_doubles(inst[1], multi[1], 1e-6));
    REQUIRE(equal_doubles(inst[2], multi[2], 1e-6));
  }
}

TEST_CASE("multi causal splitting with many identical treatments on NaN data is identical to instrumental splitting", "[multi_causal], [splitting], [NaN]") {
  auto data_vec = load_data("test/forest/resources/causal_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index({11, 11, 11, 11, 11});
  data.set_instrument_index(11);
  size_t num_features = 10;
  size_t num_treatments = data.get_num_treatments();

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto inst_splitting_rule = std::unique_ptr<SplittingRule>(new InstrumentalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty()));
  auto multi_splitting_rule = std::unique_ptr<SplittingRule>(new MultiCausalSplittingRule(
      data.get_num_rows(),
      options.get_min_node_size(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      1,
      num_treatments));

  for (auto size_node : {1000, 500, 250, 113, 25, 10, 3}) {
    std::vector<double> inst = run_splits_multi(data, size_node, options, inst_splitting_rule, relabeling_strategy, num_features);
    std::vector<double> multi = run_splits_multi(data, size_node, options, multi_splitting_rule, relabeling_strategy, num_features);

    REQUIRE(equal_doubles(inst[0], multi[0], 1e-6));
    REQUIRE(equal_doubles(inst[1], multi[1], 1e-6));
    REQUIRE(equal_doubles(inst[2], multi[2], 1e-6));
  }
}
