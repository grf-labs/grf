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

#ifndef GRF_FORESTTRAINER_H
#define GRF_FORESTTRAINER_H

#include "prediction/OptimizedPredictionStrategy.h"
#include "relabeling/RelabelingStrategy.h"
#include "splitting/factory/SplittingRuleFactory.h"

#include "tree/Tree.h"
#include "tree/TreeTrainer.h"
#include "forest/Forest.h"

#include <future>
#include <memory>
#include <set>
#include <thread>

class ForestTrainer {
public:
  ForestTrainer(std::unordered_map<size_t, size_t> observables,
                std::shared_ptr<RelabelingStrategy> relabeling_strategy,
                std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
                std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy);

  Forest train(Data* data);

  void init(uint mtry,
            uint num_trees,
            uint seed,
            uint num_threads,
            uint min_node_size,
            std::set<size_t> no_split_variables,
            std::string split_select_weights_file,
            bool sample_with_replacement,
            std::string case_weights_file,
            double sample_fraction,
            bool honesty,
            uint ci_group_size);

private:
  std::vector<std::shared_ptr<Tree>> train_ci_group(Data* data,
                                                    const Observations& observations,
                                                    RandomSampler& sampler,
                                                    double sample_fraction);

  std::vector<std::shared_ptr<Tree>> train_batch(
      uint thread_idx,
      Data* data,
      const Observations& observations);

  void set_split_select_weights(std::vector<double>& split_select_weights,
                                size_t num_independent_variables);

  size_t num_trees;
  uint mtry;
  uint min_node_size;
  uint seed;
  bool sample_with_replacement;
  double sample_fraction;

  std::set<size_t> no_split_variables;

  uint num_threads;
  std::vector<uint> thread_ranges;

  std::shared_ptr<TreeTrainer> tree_trainer;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<size_t> deterministic_vars;
  std::vector<size_t> split_select_vars;
  std::vector<double> split_select_weights;


  std::unordered_map<size_t, size_t> observables;
  std::shared_ptr<RelabelingStrategy> relabeling_strategy;
  std::shared_ptr<SplittingRuleFactory> splitting_rule_factory;
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy;

  std::vector<double> sample_weights;
  std::string sample_weights_file;
  std::string split_select_weights_file;

  uint ci_group_size;
};


#endif //GRF_FORESTTRAINER_H
