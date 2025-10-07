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

#include "splitting/factory/InstrumentalSplittingRuleFactory.h"
#include "splitting/factory/ProbabilitySplittingRuleFactory.h"
#include "splitting/factory/RegressionSplittingRuleFactory.h"
#include "splitting/factory/SurvivalSplittingRuleFactory.h"
#include "splitting/factory/CausalSurvivalSplittingRuleFactory.h"
#include "relabeling/NoopRelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "relabeling/CausalSurvivalRelabelingStrategy.h"

#include "utilities/ForestTestUtilities.h"

#include "catch.hpp"

using namespace grf;

// Splitting rule input setup to emulate one run of node zero (all data) splitting on all features
void run_one_split(const Data& data,
                   const TreeOptions& options,
                   const std::unique_ptr<SplittingRuleFactory>& splitting_rule_factory,
                   const std::unique_ptr<RelabelingStrategy>& relabeling_strategy,
                   size_t num_features,
                   size_t& split_var,
                   double& split_value) {
  std::unique_ptr<SplittingRule> splitting_rule = splitting_rule_factory->create(data.get_num_rows(), data, options);
  std::vector<size_t> possible_split_vars(num_features - 1);
  // Fill with {0, 1, 2, ..., Xj}
  std::iota(possible_split_vars.begin(), possible_split_vars.end(), 0);
  size_t node = 0;
  size_t size_node = data.get_num_rows();
  Eigen::ArrayXXd responses_by_sample(size_node, 1);
  std::vector<std::vector<size_t>> samples(1);
  for (size_t sample = 0; sample < size_node; ++sample) {
    samples[node].push_back(sample);
  }
  relabeling_strategy->relabel(samples[node], data, responses_by_sample);

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
  split_var = split_vars[node];
  split_value = split_values[node];
}

TEST_CASE("regression splitting on Xij then setting all values to the left to NaN yields the same split", "[NaN], [regression], [splitting]") {
  auto data_vec = load_data("test/forest/resources/regression_data.csv");
  Data data(data_vec);
  size_t num_features = 10;
  data.set_outcome_index(num_features);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new RegressionSplittingRuleFactory());

  size_t split_var, split_var_nan;
  double split_val, split_val_nan;
  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var, split_val);

  // Set all values to the left of the split to missing
  for(size_t row = 0; row < data.get_num_rows(); ++row) {
    double value = data.get(row, split_var);
    if (value < split_val) {
      set_data(data_vec, row, split_var, NAN);
    }
  }

  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var_nan, split_val_nan);
  REQUIRE(split_var == split_var_nan);
  REQUIRE(split_val == split_val_nan);
}

TEST_CASE("instrumental splitting on Xij then setting all values to the left to NaN yields the same split", "[NaN], [causal], [splitting]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  size_t num_features = 10;
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();
  double reduced_form_weight = 0.0;

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(reduced_form_weight));
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new InstrumentalSplittingRuleFactory());

  size_t split_var, split_var_nan;
  double split_val, split_val_nan;
  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var, split_val);

  // Set all values to the left of the split to missing
  for(size_t row = 0; row < data.get_num_rows(); ++row) {
    double value = data.get(row, split_var);
    if (value < split_val) {
      set_data(data_vec, row, split_var, NAN);
    }
  }

  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var_nan, split_val_nan);
  REQUIRE(split_var == split_var_nan);
  REQUIRE(split_val == split_val_nan);
}

TEST_CASE("probability splitting on Xij then setting all values to the left to NaN yields the same split", "[NaN], [quantile], [splitting]") {
  auto data_vec = load_data("test/forest/resources/quantile_data.csv");
  Data data(data_vec);
  size_t num_features = 10;
  data.set_outcome_index(10);
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  size_t num_classes = quantiles.size() + 1;

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new ProbabilitySplittingRuleFactory(quantiles.size() + 1));

  size_t split_var, split_var_nan;
  double split_val, split_val_nan;
  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var, split_val);

  // Set all values to the left of the split to missing
  for(size_t row = 0; row < data.get_num_rows(); ++row) {
    double value = data.get(row, split_var);
    if (value < split_val) {
      set_data(data_vec, row, split_var, NAN);
    }
  }

  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var_nan, split_val_nan);
  REQUIRE(split_var == split_var_nan);
  REQUIRE(split_val == split_val_nan);
}

TEST_CASE("survival splitting on Xij then setting all values to the left to NaN yields the same split", "[NaN], [survival], [splitting]") {
  auto data_vec = load_data("test/forest/resources/survival_data_MIA.csv");
  Data data(data_vec);
  size_t num_features = 5;
  data.set_outcome_index(5);
  data.set_censor_index(6);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  bool fast_logrank = false;
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new SurvivalSplittingRuleFactory(fast_logrank));

  size_t split_var, split_var_nan;
  double split_val, split_val_nan;
  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var, split_val);

  // Set all values to the left of the split to missing
  for(size_t row = 0; row < data.get_num_rows(); ++row) {
    double value = data.get(row, split_var);
    if (value < split_val) {
      set_data(data_vec, row, split_var, NAN);
    }
  }

  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var_nan, split_val_nan);
  REQUIRE(split_var == split_var_nan);
  REQUIRE(split_val == split_val_nan);
}

TEST_CASE("accelerated survival splitting on Xij then setting all values to the left to NaN yields the same split", "[NaN], [survival], [splitting]") {
  auto data_vec = load_data("test/forest/resources/survival_data_MIA.csv");
  Data data(data_vec);
  size_t num_features = 5;
  data.set_outcome_index(5);
  data.set_censor_index(6);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  bool fast_logrank = true;
  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new SurvivalSplittingRuleFactory(fast_logrank));

  size_t split_var, split_var_nan;
  double split_val, split_val_nan;
  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var, split_val);

  // Set all values to the left of the split to missing
  for(size_t row = 0; row < data.get_num_rows(); ++row) {
    double value = data.get(row, split_var);
    if (value < split_val) {
      set_data(data_vec, row, split_var, NAN);
    }
  }

  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var_nan, split_val_nan);
  REQUIRE(split_var == split_var_nan);
  REQUIRE(split_val == split_val_nan);
}

TEST_CASE("causal survival splitting on Xij then setting all values to the left to NaN yields the same split", "[NaN], [causal survival], [splitting]") {
  auto data_vec = load_data("test/forest/resources/causal_survival_data.csv");
  Data data(data_vec);
  size_t num_features = 5;
  data.set_treatment_index(5);
  data.set_instrument_index(5);
  data.set_censor_index(6);
  data.set_causal_survival_numerator_index(7);
  data.set_causal_survival_denominator_index(8);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new CausalSurvivalRelabelingStrategy());
  auto splitting_rule_factory = std::unique_ptr<SplittingRuleFactory>(new CausalSurvivalSplittingRuleFactory());

  size_t split_var, split_var_nan;
  double split_val, split_val_nan;
  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var, split_val);

  // Set all values to the left of the split to missing
  for(size_t row = 0; row < data.get_num_rows(); ++row) {
    double value = data.get(row, split_var);
    if (value < split_val) {
      set_data(data_vec, row, split_var, NAN);
    }
  }

  run_one_split(data, options, splitting_rule_factory, relabeling_strategy, num_features, split_var_nan, split_val_nan);
  REQUIRE(split_var == split_var_nan);
  REQUIRE(split_val == split_val_nan);
}
