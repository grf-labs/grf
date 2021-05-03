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

#include "commons/utility.h"
#include "ForestTrainer.h"
#include "random/random.hpp"


namespace grf {

ForestTrainer::ForestTrainer(std::unique_ptr<RelabelingStrategy> relabeling_strategy,
                             std::unique_ptr<SplittingRuleFactory> splitting_rule_factory,
                             std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy) :
    tree_trainer(std::move(relabeling_strategy),
                 std::move(splitting_rule_factory),
                 std::move(prediction_strategy)) {}

Forest ForestTrainer::train(const Data& data, const ForestOptions& options) const {
  std::vector<std::unique_ptr<Tree>> trees = train_trees(data, options);

  size_t num_variables = data.get_num_cols() - data.get_disallowed_split_variables().size();
  size_t ci_group_size = options.get_ci_group_size();
  return Forest(trees, num_variables, ci_group_size);
}

std::vector<std::unique_ptr<Tree>> ForestTrainer::train_trees(const Data& data,
                                                              const ForestOptions& options) const {
  size_t num_samples = data.get_num_rows();
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

  std::vector<std::future<std::vector<std::unique_ptr<Tree>>>> futures;
  futures.reserve(thread_ranges.size());

  std::vector<std::unique_ptr<Tree>> trees;
  trees.reserve(num_trees);

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_trees_batch = thread_ranges[i + 1] - start_index;

    futures.push_back(std::async(std::launch::async,
                                 &ForestTrainer::train_batch,
                                 this,
                                 start_index,
                                 num_trees_batch,
                                 std::ref(data),
                                 options));
  }

  for (auto& future : futures) {
    std::vector<std::unique_ptr<Tree>> thread_trees = future.get();
    trees.insert(trees.end(),
                 std::make_move_iterator(thread_trees.begin()),
                 std::make_move_iterator(thread_trees.end()));
  }

  return trees;
}

std::vector<std::unique_ptr<Tree>> ForestTrainer::train_batch(
    size_t start,
    size_t num_trees,
    const Data& data,
    const ForestOptions& options) const {
  size_t ci_group_size = options.get_ci_group_size();

  std::mt19937_64 random_number_generator(options.get_random_seed() + start);
  nonstd::uniform_int_distribution<uint> udist;
  std::vector<std::unique_ptr<Tree>> trees;
  trees.reserve(num_trees * ci_group_size);

  for (size_t i = 0; i < num_trees; i++) {
    uint tree_seed = udist(random_number_generator);
    RandomSampler sampler(tree_seed, options.get_sampling_options());

    if (ci_group_size == 1) {
      std::unique_ptr<Tree> tree = train_tree(data, sampler, options);
      trees.push_back(std::move(tree));
    } else {
      std::vector<std::unique_ptr<Tree>> group = train_ci_group(data, sampler, options);
      trees.insert(trees.end(),
          std::make_move_iterator(group.begin()),
          std::make_move_iterator(group.end()));
    }
  }
  return trees;
}
std::unique_ptr<Tree> ForestTrainer::train_tree(const Data& data,
                                                RandomSampler& sampler,
                                                const ForestOptions& options) const {
  std::vector<size_t> clusters;
  sampler.sample_clusters(data.get_num_rows(), options.get_sample_fraction(), clusters);
  return tree_trainer.train(data, sampler, clusters, options.get_tree_options());
}

std::vector<std::unique_ptr<Tree>> ForestTrainer::train_ci_group(const Data& data,
                                                                 RandomSampler& sampler,
                                                                 const ForestOptions& options) const {
  std::vector<std::unique_ptr<Tree>> trees;

  std::vector<size_t> clusters;
  sampler.sample_clusters(data.get_num_rows(), 0.5, clusters);

  double sample_fraction = options.get_sample_fraction();
  for (size_t i = 0; i < options.get_ci_group_size(); ++i) {
    std::vector<size_t> cluster_subsample;
    sampler.subsample(clusters, sample_fraction * 2, cluster_subsample);

    std::unique_ptr<Tree> tree = tree_trainer.train(data, sampler, cluster_subsample, options.get_tree_options());
    trees.push_back(std::move(tree));
  }
  return trees;
}

} // namespace grf
