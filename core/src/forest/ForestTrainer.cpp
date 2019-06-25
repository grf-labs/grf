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

#include <algorithm>
#include <ctime>
#include <future>
#include <stdexcept>
#include <string>

#include "commons/utility.h"
#include "ForestTrainer.h"

ForestTrainer::ForestTrainer(std::shared_ptr<RelabelingStrategy> relabeling_strategy,
                             std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
                             std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy) :
    tree_trainer(relabeling_strategy,
                 splitting_rule_factory,
                 prediction_strategy) {}

const Forest ForestTrainer::train(const Data* data,
                                  const ForestOptions& options) const {
  size_t num_samples = data->get_num_rows();
  uint num_trees = options.get_num_trees();

  // Ensure that the sample fraction is not too small and honesty fraction is not too extreme.
  const TreeOptions& tree_options = options.get_tree_options();
  bool honesty = tree_options.get_honesty();
  double honesty_fraction = tree_options.get_honesty_fraction();
  if ((size_t) num_samples * options.get_sample_fraction() < 1) {
    throw std::runtime_error("The sample fraction is too small, as no observations will be sampled.");
  } else if (honesty && ((size_t) num_samples * options.get_sample_fraction() * honesty_fraction < 1
             || (size_t) num_samples * options.get_sample_fraction() * (1-honesty_fraction) < 1)) {
    throw std::runtime_error("The honesty fraction is too close to 1 or 0, as no observations will be sampled.");
  }

  uint num_groups = (uint) num_trees / options.get_ci_group_size();

  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, num_groups - 1, options.get_num_threads());

  std::vector<std::future<std::vector<std::shared_ptr<Tree>>>> futures;
  futures.reserve(thread_ranges.size());

  std::vector<std::shared_ptr<Tree>> trees;
  trees.reserve(num_trees);

  std::mt19937_64 thread_seed_engine(options.get_random_seed());

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_trees_batch = thread_ranges[i + 1] - start_index;
    size_t batch_seed = thread_seed_engine();
    futures.push_back(std::async(std::launch::async,
                                 &ForestTrainer::train_batch,
                                 this,
                                 batch_seed,
                                 num_trees_batch,
                                 data,
                                 options));
  }

  for (auto& future : futures) {
    std::vector<std::shared_ptr<Tree>> thread_trees = future.get();
    trees.insert(trees.end(), thread_trees.begin(), thread_trees.end());
  }
  return Forest::create(trees, options, data);
}

std::vector<std::shared_ptr<Tree>> ForestTrainer::train_batch(
    size_t batch_seed,
    size_t num_trees,
    const Data* data,
    const ForestOptions& options) const {
  size_t ci_group_size = options.get_ci_group_size();

  std::vector<std::shared_ptr<Tree>> trees;

  if (ci_group_size == 1) {
    trees.reserve(num_trees);
  } else {
    trees.reserve(num_trees * ci_group_size);
  }

  for (size_t i = 0; i < num_trees; i++) {
    uint tree_seed = batch_seed + i;
    RandomSampler sampler(tree_seed, options.get_sampling_options());

    if (ci_group_size == 1) {
      std::shared_ptr<Tree> tree = train_tree(data, sampler, options);
      trees.push_back(tree);
    } else {
      std::vector<std::shared_ptr<Tree>> group = train_ci_group(data, sampler, options);
      trees.insert(trees.end(), group.begin(), group.end());
    }
  }
  return trees;
}

std::shared_ptr<Tree> ForestTrainer::train_tree(const Data* data,
                                                RandomSampler& sampler,
                                                const ForestOptions& options) const {
  std::vector<size_t> clusters;
  sampler.sample_clusters(data->get_num_rows(), options.get_sample_fraction(), clusters);
  return tree_trainer.train(data, sampler, clusters, options.get_tree_options());
}

std::vector<std::shared_ptr<Tree>> ForestTrainer::train_ci_group(const Data* data,
                                                                 RandomSampler& sampler,
                                                                 const ForestOptions& options) const {
  std::vector<std::shared_ptr<Tree>> trees;

  std::vector<size_t> clusters;
  sampler.sample_clusters(data->get_num_rows(), 0.5, clusters);

  double sample_fraction = options.get_sample_fraction();
  for (size_t i = 0; i < options.get_ci_group_size(); ++i) {
    std::vector<size_t> cluster_subsample;
    sampler.subsample(clusters, sample_fraction * 2, cluster_subsample);

    std::shared_ptr<Tree> tree = tree_trainer.train(data,
        sampler, cluster_subsample, options.get_tree_options());
    trees.push_back(tree);
  }
  return trees;
}
