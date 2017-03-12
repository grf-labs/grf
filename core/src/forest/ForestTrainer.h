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

#ifndef GRADIENTFOREST_FORESTTRAINER_H
#define GRADIENTFOREST_FORESTTRAINER_H

#include "RelabelingStrategy.h"
#include "SplittingRuleFactory.h"
#include "PredictionStrategy.h"

#include "Tree.h"
#include "TreeTrainer.h"
#include "Forest.h"

#include <memory>
#include <thread>
#include <future>

class ForestTrainer {
public:
  ForestTrainer(std::unordered_map<size_t, size_t> observables,
                std::shared_ptr<RelabelingStrategy> relabeling_strategy,
                std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
                std::shared_ptr<PredictionStrategy> prediction_strategy);

  Forest train(Data* data);

  void init(uint mtry,
            uint num_trees, std::ostream *verbose_out, uint seed, uint num_threads,
            std::string load_forest_filename, uint min_node_size,
            std::vector<size_t> no_split_variables,
            std::string split_select_weights_file, std::vector<std::string>& always_split_variable_names,
            bool sample_with_replacement,
            bool memory_saving_splitting,
            std::string case_weights_file,
            double sample_fraction,
            bool honesty,
            uint ci_group_size);

private:
  std::vector<std::shared_ptr<Tree>> train_ci_group(Data* data,
                                                    const Observations& observations,
                                                    BootstrapSampler& bootstrap_sampler,
                                                    double sample_fraction);

  std::vector<std::shared_ptr<Tree>> train_batch(
      uint thread_idx,
      Data* data,
      const Observations& observations);

  void set_split_select_weights(std::vector<double>& split_select_weights,
                                size_t num_independent_variables);
  void set_always_split_variables(Data* data,
                                  std::vector<std::string>& always_split_variable_names,
                                  size_t num_independent_variables);

  std::ostream* verbose_out;

  size_t num_trees;
  uint mtry;
  uint min_node_size;
  uint seed;
  bool prediction_mode;
  bool sample_with_replacement;
  bool memory_saving_splitting;
  double sample_fraction;

  std::vector<size_t> no_split_variables;

  // Multithreading
  uint num_threads;
  std::vector<uint> thread_ranges;

  TreeTrainer* tree_trainer;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<std::string> always_split_variable_names;
  std::vector<size_t> deterministic_varIDs;
  std::vector<size_t> split_select_varIDs;
  std::vector<double> split_select_weights;

  // Bootstrap weights
  std::vector<double> case_weights;

  std::unordered_map<size_t, size_t> observables;
  std::shared_ptr<RelabelingStrategy> relabeling_strategy;
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory;
  std::shared_ptr<PredictionStrategy> prediction_strategy;

  std::string split_select_weights_file;
  std::string case_weights_file;

  uint ci_group_size;
};


#endif //GRADIENTFOREST_FORESTTRAINER_H
