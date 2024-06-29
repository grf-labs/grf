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
#include "splitting/RegressionSplittingRule.h"
#include "splitting/MultiRegressionSplittingRule.h"
#include "relabeling/NoopRelabelingStrategy.h"

#include "commons/utility.h"
#include "utilities/ForestTestUtilities.h"
#include "utilities/FileTestUtilities.h"

#include "catch.hpp"

using namespace grf;

// Splitting rule input setup to emulate one run of node zero (all data) splitting on all features
// returning a vector containing the best split variable, best split value, and missing direction.
std::vector<double> run_splits(const Data& data,
                               const TreeOptions& options,
                               const std::unique_ptr<SplittingRule>& splitting_rule,
                               const std::unique_ptr<RelabelingStrategy>& relabeling_strategy,
                               size_t num_features) {
  size_t node = 0;
  size_t size_node = data.get_num_rows();
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

TEST_CASE("multi regression splitting with one outcome is identical to regression splitting", "[multi_regression], [splitting]") {
  auto data_vec = load_data("test/forest/resources/regression_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  size_t num_features = 10;

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto reg_splitting_rule = std::unique_ptr<SplittingRule>(new RegressionSplittingRule(
      data.get_num_rows(),
      options.get_alpha(),
      options.get_imbalance_penalty()));
  auto multi_reg_splitting_rule = std::unique_ptr<SplittingRule>(new MultiRegressionSplittingRule(
      data.get_num_rows(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      1));


  std::vector<double> reg = run_splits(data, options, reg_splitting_rule, relabeling_strategy, num_features);
  std::vector<double> multi_reg = run_splits(data, options, multi_reg_splitting_rule, relabeling_strategy, num_features);

  REQUIRE(equal_doubles(reg[0], multi_reg[0], 1e-6));
  REQUIRE(equal_doubles(reg[1], multi_reg[1], 1e-6));
  REQUIRE(equal_doubles(reg[2], multi_reg[2], 1e-6));
}

TEST_CASE("multi regression splitting with one outcome on NaN data is identical to regression splitting", "[multi_regression], [splitting], [NaN]") {
  auto data_vec = load_data("test/forest/resources/regression_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(5);
  size_t num_features = 5;

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto reg_splitting_rule = std::unique_ptr<SplittingRule>(new RegressionSplittingRule(
      data.get_num_rows(),
      options.get_alpha(),
      options.get_imbalance_penalty()));
  auto multi_reg_splitting_rule = std::unique_ptr<SplittingRule>(new MultiRegressionSplittingRule(
      data.get_num_rows(),
      options.get_alpha(),
      options.get_imbalance_penalty(),
      1));


  std::vector<double> reg = run_splits(data, options, reg_splitting_rule, relabeling_strategy, num_features);
  std::vector<double> multi_reg = run_splits(data, options, multi_reg_splitting_rule, relabeling_strategy, num_features);

  REQUIRE(equal_doubles(reg[0], multi_reg[0], 1e-6));
  REQUIRE(equal_doubles(reg[1], multi_reg[1], 1e-6));
  REQUIRE(equal_doubles(reg[2], multi_reg[2], 1e-6));
}
