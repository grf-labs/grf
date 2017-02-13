/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <memory>

#include "Data.h"
#include "TreeTrainer.h"
#include "TreePredictor.h"

TreeTrainer::TreeTrainer(std::shared_ptr<RelabelingStrategy> relabeling_strategy,
                         std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
                         std::shared_ptr<PredictionStrategy> prediction_strategy,
                         const TreeOptions& options) :
    relabeling_strategy(relabeling_strategy),
    splitting_rule_factory(splitting_rule_factory),
    prediction_strategy(prediction_strategy),
    options(options) {}

std::shared_ptr<Tree> TreeTrainer::train(Data* data,
                                         const Observations& observations,
                                         BootstrapSampler& bootstrap_sampler,
                                         const std::vector<size_t>& sampleIDs) {
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> nodes;
  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;

  child_nodeIDs.push_back(std::vector<size_t>());
  child_nodeIDs.push_back(std::vector<size_t>());
  create_empty_node(child_nodeIDs, nodes, split_varIDs, split_values);

  std::vector<size_t> leaf_sampleIDs;

  if (options.get_honesty()) {
    bootstrap_sampler.subsample(sampleIDs, 0.5, nodes[0], leaf_sampleIDs);
  } else {
    nodes[0] = sampleIDs;
  }

  std::shared_ptr<SplittingRule> splitting_rule = splitting_rule_factory->create();

  size_t num_open_nodes = 1;
  size_t i = 0;
  while (num_open_nodes > 0) {
    bool is_terminal_node = split_node(i,
                                       splitting_rule,
                                       bootstrap_sampler,
                                       data,
                                       observations,
                                       child_nodeIDs,
                                       nodes,
                                       split_varIDs,
                                       split_values,
                                       options.get_split_select_weights());
    if (is_terminal_node) {
      --num_open_nodes;
    } else {
      nodes[i].clear();
      ++num_open_nodes;
    }
    ++i;
  }

  auto tree = std::shared_ptr<Tree>(new Tree(child_nodeIDs,
      nodes,
      split_varIDs,
      split_values,
      std::vector<size_t>(),
      PredictionValues()));

  if (!leaf_sampleIDs.empty()) {
    repopulate_terminal_nodeIDs(tree, data, leaf_sampleIDs);
  }

  PredictionValues prediction_values = prediction_strategy->precompute_prediction_values(
      tree->get_leaf_nodeIDs(), observations);
  tree->set_prediction_values(prediction_values);

  return tree;
}

void TreeTrainer::repopulate_terminal_nodeIDs(std::shared_ptr<Tree> tree,
                                              Data* data,
                                              const std::vector<size_t>& leaf_sampleIDs) {
  size_t num_nodes = tree->get_leaf_nodeIDs().size();
  std::vector<std::vector<size_t>> new_terminal_nodeIDs(num_nodes);

  TreePredictor tree_predictor;
  std::vector<size_t> terminal_nodeIDs = tree_predictor.get_terminal_nodeIDs(
      tree, data, leaf_sampleIDs);

  for (auto& sampleID : leaf_sampleIDs) {
    size_t terminal_nodeID = terminal_nodeIDs.at(sampleID);
    new_terminal_nodeIDs.at(terminal_nodeID).push_back(sampleID);
  }
  tree->set_leaf_nodeIDs(new_terminal_nodeIDs);
}

void TreeTrainer::create_split_variable_subset(std::vector<size_t>& result,
                                               BootstrapSampler &bootstrap_sampler,
                                               Data *data,
                                               const std::vector<double>& split_select_weights) {

  // Always use deterministic variables
  std::vector<size_t> deterministic_varIDs = options.get_deterministic_varIDs();
  std::copy(deterministic_varIDs.begin(), deterministic_varIDs.end(), std::inserter(result, result.end()));

  // Randomly add non-deterministic variables (according to weights if needed)
  uint mtry = options.get_mtry();
  if (split_select_weights.empty()) {
    bootstrap_sampler.draw_without_replacement_skip(result,
                                                    data->get_num_cols(),
                                                    options.get_no_split_variables(),
                                                    mtry);
  } else {
    size_t num_draws = mtry - result.size();
    bootstrap_sampler.draw_without_replacement_weighted(result,
                                                        options.get_split_select_varIDs(),
                                                        num_draws,
                                                        split_select_weights);
  }
}

bool TreeTrainer::split_node(size_t nodeID,
                             std::shared_ptr<SplittingRule> splitting_rule,
                             BootstrapSampler& bootstrap_sampler,
                             Data* data,
                             const Observations& observations,
                             std::vector<std::vector<size_t>>& child_nodeIDs,
                             std::vector<std::vector<size_t>>& sampleIDs,
                             std::vector<size_t>& split_varIDs,
                             std::vector<double>& split_values,
                             const std::vector<double>& split_select_weights) {

// Select random subset of variables to possibly split at
  std::vector<size_t> possible_split_varIDs;
  create_split_variable_subset(possible_split_varIDs, bootstrap_sampler, data, split_select_weights);

// Call subclass method, sets split_varIDs and split_values
  bool stop = split_node_internal(nodeID,
                                  splitting_rule,
                                  observations,
                                  possible_split_varIDs,
                                  sampleIDs,
                                  split_varIDs,
                                  split_values);
  if (stop) {
    // Terminal node
    return true;
  }

  size_t split_varID = split_varIDs[nodeID];
  double split_value = split_values[nodeID];

// Create child nodes
  size_t left_child_nodeID = sampleIDs.size();
  child_nodeIDs[0][nodeID] = left_child_nodeID;
  create_empty_node(child_nodeIDs, sampleIDs, split_varIDs, split_values);

  size_t right_child_nodeID = sampleIDs.size();
  child_nodeIDs[1][nodeID] = right_child_nodeID;
  create_empty_node(child_nodeIDs, sampleIDs, split_varIDs, split_values);

  // For each sample in node, assign to left or right child
  // Ordered: left is <= splitval and right is > splitval
  for (auto& sampleID : sampleIDs[nodeID]) {
    if (data->get(sampleID, split_varID) <= split_value) {
      sampleIDs[left_child_nodeID].push_back(sampleID);
    } else {
      sampleIDs[right_child_nodeID].push_back(sampleID);
    }
  }

  // No terminal node
  return false;
}

bool TreeTrainer::split_node_internal(size_t nodeID,
                                      std::shared_ptr<SplittingRule> splitting_rule,
                                      const Observations& observations,
                                      const std::vector<size_t>& possible_split_varIDs,
                                      std::vector<std::vector<size_t>>& sampleIDs,
                                      std::vector<size_t>& split_varIDs,
                                      std::vector<double>& split_values) {
  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= options.get_min_node_size()) {
    split_values[nodeID] = -1.0;
    return true;
  }

  // Check if node is pure and set split_value to estimate and stop if pure
  bool pure = true;
  double pure_value = 0;
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    size_t sampleID = sampleIDs[nodeID][i];
    double value = observations.get(Observations::OUTCOME, sampleID);
    if (i != 0 && value != pure_value) {
      pure = false;
      break;
    }
    pure_value = value;
  }

  if (pure) {
    split_values[nodeID] = -1.0;
    return true;
  }

  std::unordered_map<size_t, double> responses_by_sampleID = relabeling_strategy->relabel_outcomes(
      observations, sampleIDs[nodeID]);

  bool stop = responses_by_sampleID.empty() ||
              splitting_rule->find_best_split(nodeID,
                                              possible_split_varIDs,
                                              responses_by_sampleID,
                                              sampleIDs,
                                              split_varIDs,
                                              split_values);

  if (stop) {
    split_values[nodeID] = -1.0;
    return true;
  }
  return false;
}

void TreeTrainer::create_empty_node(std::vector<std::vector<size_t>>& child_nodeIDs,
                                    std::vector<std::vector<size_t>>& sampleIDs,
                                    std::vector<size_t>& split_varIDs,
                                    std::vector<double>& split_values) {
  split_varIDs.push_back(0);
  split_values.push_back(0);
  child_nodeIDs[0].push_back(0);
  child_nodeIDs[1].push_back(0);
  sampleIDs.push_back(std::vector<size_t>());
}
