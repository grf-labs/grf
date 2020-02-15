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

#include "catch.hpp"

#define private public // to access and test the private members `find_best_split_value_small_q`
                       // and `find_best_split_value_big_q` of a splitting rule.

#include "commons/utility.h"
#include "relabeling/NoopRelabelingStrategy.h"
#include "relabeling/InstrumentalRelabelingStrategy.h"
#include "relabeling/QuantileRelabelingStrategy.h"
#include "splitting/RegressionSplittingRule.h"
#include "splitting/InstrumentalSplittingRule.h"
#include "splitting/ProbabilitySplittingRule.h"
#include "utilities/ForestTestUtilities.h"

using namespace grf;

TEST_CASE("regression splitting small q / big Q are identical", "[regression, splitting]") {
  std::unique_ptr<Data> data = load_data("test/forest/resources/regression_data.csv");
  size_t weight_index = 9;
  data->set_weight_index(weight_index);
  data->set_outcome_index(10);

  // Use covariate in data column 9 as dummy sample weights
  bool error;
  for(size_t r = 0; r < data->get_num_rows(); r++) {
    double value = data->get(r, weight_index);
    double weight = value < 0 ? -value : value;
    data->set(weight_index, r, weight, error);
  }

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new NoopRelabelingStrategy());
  auto splitting_rule = std::unique_ptr<RegressionSplittingRule>(new RegressionSplittingRule(
    data->get_max_num_unique_values(),
    options.get_alpha(),
    options.get_imbalance_penalty()));

  // Splitting rule input setup to emulate one run of node zero (all data) splitting on all features
  std::vector<size_t> possible_split_vars(weight_index - 1);
  std::iota(possible_split_vars.begin(), possible_split_vars.end(), 0); // {0, 1, 2, ..., Xj}

  size_t node = 0;
  double weight_sum_node = 0;
  double sum_node = 0;
  size_t size_node = data->get_num_rows();
  size_t min_child_size = 1;
  std::vector<double> responses_by_sample (size_node);
  std::vector<std::vector<size_t>> samples (1);

  for (size_t sample = 0; sample < size_node; ++sample) {
    samples[node].push_back(sample);
  }

  relabeling_strategy->relabel(samples[node], *data, responses_by_sample);

  for (auto& sample : samples[node]) {
    double sample_weight = data->get_weight(sample);
    weight_sum_node += sample_weight;
    sum_node += sample_weight * responses_by_sample[sample];
  }

  double best_value = 0;
  size_t best_var = 0;
  double best_decrease = 0;
  bool best_send_missing_left = true;

  // Assert all splits are the same
  for (auto& var : possible_split_vars) {
    splitting_rule->find_best_split_value_small_q(*data, node, var, weight_sum_node, sum_node, size_node, min_child_size,
                                                  best_value, best_var, best_decrease, best_send_missing_left, responses_by_sample, samples);
    size_t best_var_q = best_var;
    double best_value_q = best_value;
    double best_decrease_q = best_decrease;

    best_value = 0;
    best_var = 0;
    best_decrease = 0;

    splitting_rule->find_best_split_value_large_q(*data, node, var, weight_sum_node, sum_node, size_node, min_child_size,
                                                  best_value, best_var, best_decrease, responses_by_sample, samples);
    REQUIRE(equal_doubles(best_value_q, best_value, 1e-9));
    REQUIRE(best_var_q == best_var);
    REQUIRE(equal_doubles(best_decrease_q, best_decrease, 1e-9));
    best_value = 0;
    best_var = 0;
    best_decrease = 0;
  }
}

TEST_CASE("instrumental splitting small q / big Q are identical", "[causal, splitting]") {
  std::unique_ptr<Data> data = load_data("test/forest/resources/causal_data.csv");
  size_t weight_index = 9;
  data->set_weight_index(weight_index);
  data->set_outcome_index(10);
  data->set_treatment_index(11);
  data->set_instrument_index(11);

  // Use covariate in data column 9 as dummy sample weights
  bool error;
  for(size_t r = 0; r < data->get_num_rows(); r++) {
    double value = data->get(r, weight_index);
    double weight = value < 0 ? -value : value;
    data->set(weight_index, r, weight, error);
  }

  double reduced_form_weight = 0.0;

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new InstrumentalRelabelingStrategy(reduced_form_weight));
  auto splitting_rule = std::unique_ptr<InstrumentalSplittingRule>(new InstrumentalSplittingRule(
    data->get_max_num_unique_values(),
    options.get_min_node_size(),
    options.get_alpha(),
    options.get_imbalance_penalty()));

  // Splitting rule input setup to emulate one run of node zero (all data) splitting on all features
  std::vector<size_t> possible_split_vars(weight_index - 1);
  std::iota(possible_split_vars.begin(), possible_split_vars.end(), 0); // {0, 1, 2, ..., Xj}

  size_t node = 0;
  double weight_sum_node = 0.0;
  double sum_node = 0.0;
  double sum_node_z = 0.0;
  double sum_node_z_squared = 0.0;
  size_t num_samples = data->get_num_rows();
  std::vector<double> responses_by_sample (num_samples);
  std::vector<std::vector<size_t>> samples (1);

  for (size_t sample = 0; sample < num_samples; ++sample) {
    samples[node].push_back(sample);
  }

  relabeling_strategy->relabel(samples[node], *data, responses_by_sample);

  for (auto& sample : samples[node]) {
    double sample_weight = data->get_weight(sample);
    weight_sum_node += sample_weight;
    sum_node += sample_weight * responses_by_sample[sample];

    double z = data->get_instrument(sample);
    sum_node_z += sample_weight * z;
    sum_node_z_squared += sample_weight * z * z;
  }

  double size_node = sum_node_z_squared - sum_node_z * sum_node_z / weight_sum_node;
  double min_child_size = size_node * options.get_alpha();

  double mean_z_node = sum_node_z / weight_sum_node;
  size_t num_node_small_z = 0;
  for (auto& sample : samples[node]) {
    double z = data->get_instrument(sample);
    if (z < mean_z_node) {
      num_node_small_z++;
    }
  }

  size_t best_var = 0;
  double best_value = 0;
  double best_decrease = 0.0;
  bool best_send_missing_left = true;

  // Assert all splits are the same
  for (auto& var : possible_split_vars) {
    splitting_rule->find_best_split_value_small_q(*data, node, var, num_samples, weight_sum_node, sum_node, mean_z_node, num_node_small_z,
                                                  sum_node_z, sum_node_z_squared, min_child_size, best_value,
                                                  best_var, best_decrease, best_send_missing_left, responses_by_sample, samples);
    size_t best_var_q = best_var;
    double best_value_q = best_value;
    double best_decrease_q = best_decrease;

    best_value = 0;
    best_var = 0;
    best_decrease = 0;

    splitting_rule->find_best_split_value_large_q(*data, node, var, num_samples, weight_sum_node, sum_node, mean_z_node, num_node_small_z,
                                                  sum_node_z, sum_node_z_squared, min_child_size, best_value,
                                                  best_var, best_decrease, responses_by_sample, samples);
    REQUIRE(equal_doubles(best_value_q, best_value, 1e-9));
    REQUIRE(best_var_q == best_var);
    REQUIRE(equal_doubles(best_decrease_q, best_decrease, 1e-9));
    best_value = 0;
    best_var = 0;
    best_decrease = 0;
  }
}

TEST_CASE("quantile splitting small q / big Q are identical", "[quantile, splitting]") {
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  size_t num_classes = quantiles.size() + 1;
  std::unique_ptr<Data> data = load_data("test/forest/resources/quantile_data.csv");
  size_t outcome_index = 10;
  data->set_outcome_index(outcome_index);

  TreeOptions options = ForestTestUtilities::default_options().get_tree_options();

  std::unique_ptr<RelabelingStrategy> relabeling_strategy(new QuantileRelabelingStrategy(quantiles));
  auto splitting_rule = std::unique_ptr<ProbabilitySplittingRule>(new ProbabilitySplittingRule(
    data->get_max_num_unique_values(),
    num_classes,
    options.get_alpha(),
    options.get_imbalance_penalty()));

  // Splitting rule input setup to emulate one run of node zero (all data) splitting on all features
  std::vector<size_t> possible_split_vars(outcome_index - 1);
  std::iota(possible_split_vars.begin(), possible_split_vars.end(), 0); // {0, 1, 2, ..., Xj}

  size_t node = 0;
  size_t size_node = data->get_num_rows();
  size_t min_child_size = 1;
  std::vector<double> responses_by_sample (size_node);
  std::vector<std::vector<size_t>> samples (1);

  for (size_t sample = 0; sample < size_node; ++sample) {
    samples[node].push_back(sample);
  }

  relabeling_strategy->relabel(samples[node], *data, responses_by_sample);

  size_t* class_counts = new size_t[num_classes]();
  for (size_t i = 0; i < size_node; ++i) {
    size_t sample = samples[node][i];
    uint sample_class = (uint) std::round(responses_by_sample[sample]);
    ++class_counts[sample_class];
  }

  double best_value = 0;
  size_t best_var = 0;
  double best_decrease = 0;
  bool best_send_missing_left = true;

  // Assert all splits are the same
  for (auto& var : possible_split_vars) {
    splitting_rule->find_best_split_value_small_q(*data, node, var, num_classes, class_counts, size_node, min_child_size,
                                                  best_value, best_var, best_decrease, best_send_missing_left, responses_by_sample, samples);

    size_t best_var_q = best_var;
    double best_value_q = best_value;
    double best_decrease_q = best_decrease;

    best_value = 0;
    best_var = 0;
    best_decrease = 0;

    splitting_rule->find_best_split_value_large_q(*data, node, var, num_classes, class_counts, size_node, min_child_size,
                                                  best_value, best_var, best_decrease, responses_by_sample, samples);
    REQUIRE(equal_doubles(best_value_q, best_value, 1e-9));
    REQUIRE(best_var_q == best_var);
    REQUIRE(equal_doubles(best_decrease_q, best_decrease, 1e-9));
    best_value = 0;
    best_var = 0;
    best_decrease = 0;
  }
  delete[] class_counts;
}
