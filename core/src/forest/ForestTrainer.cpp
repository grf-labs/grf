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

#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <thread>
#include <future>
#include "commons/utility.h"
#include "forest/ForestTrainer.h"
#include "splitting/factory/SplittingRuleFactory.h"

ForestTrainer::ForestTrainer(std::unordered_map<size_t, size_t> observables,
                             std::shared_ptr<RelabelingStrategy> relabeling_strategy,
                             std::shared_ptr<SplittingRuleFactory> splitting_rule_factory,
                             std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy) :
    num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), seed(0),
    sample_with_replacement(true), sample_fraction(1), num_threads(DEFAULT_NUM_THREADS),
    observables(observables), relabeling_strategy(relabeling_strategy),
    splitting_rule_factory(splitting_rule_factory),
    prediction_strategy(prediction_strategy) {}

void ForestTrainer::init(uint mtry,
                         uint num_trees,
                         uint seed,
                         uint num_threads,
                         uint min_node_size,
                         std::set<size_t> no_split_variables,
                         std::string split_select_weights_file,
                         bool sample_with_replacement,
                         std::string sample_weights_file,
                         double sample_fraction,
                         bool honesty,
                         uint ci_group_size) {
  this->split_select_weights_file = split_select_weights_file;
  this->sample_weights_file = sample_weights_file;

  // Set number of threads
  if (num_threads == DEFAULT_NUM_THREADS) {
    this->num_threads = std::thread::hardware_concurrency();
  } else {
    this->num_threads = num_threads;
  }

  // If necessary, round the number of trees up to a multiple of
  // the confidence interval group size.
  this->num_trees = num_trees + (num_trees % ci_group_size);

  this->mtry = mtry;
  this->min_node_size = min_node_size;
  this->sample_with_replacement = sample_with_replacement;

  if (ci_group_size > 1 && sample_fraction > 0.5) {
    throw std::runtime_error("When confidence intervals are enabled, the"
        " sampling fraction must be less than 0.5.");
  }
  this->sample_fraction = sample_fraction;

  this->no_split_variables = no_split_variables;
  for (auto it : observables) {
      this->no_split_variables.insert(it.second);
  }

  if (seed != 0) {
    this->seed = seed;
  } else {
    std::random_device random_device;
    this->seed = random_device();
  }

  this->ci_group_size = ci_group_size;

  TreeOptions tree_options(mtry,
      min_node_size,
      split_select_weights,
      split_select_vars,
      deterministic_vars,
      this->no_split_variables,
      honesty);
  tree_trainer = std::shared_ptr<TreeTrainer>(new TreeTrainer(
      relabeling_strategy,
      splitting_rule_factory,
      prediction_strategy,
      tree_options));
}

Forest ForestTrainer::train(Data* data) {
  size_t num_samples = data->get_num_rows();
  size_t num_variables = data->get_num_cols();
  size_t num_independent_variables = num_variables - no_split_variables.size();

  if (!split_select_weights_file.empty()) {
    std::vector<double> split_select_weights;
    split_select_weights.resize(1);
    read_vector_from_file(split_select_weights, split_select_weights_file);
    if (split_select_weights.size() != num_variables - 1) {
      throw std::runtime_error("Number of split select weights is not equal to number of independent variables.");
    }
    set_split_select_weights(split_select_weights, num_independent_variables);
  }

  // Load case weights from file
  if (!sample_weights_file.empty()) {
    read_vector_from_file(sample_weights, sample_weights_file);
    if (sample_weights.size() != num_samples - 1) {
      throw std::runtime_error("Number of case weights is not equal to number of samples.");
    }
  }

  if (mtry == 0) {
    mtry = std::ceil(num_independent_variables / 3.0);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }

  // Check if any observations samples
  if ((size_t) num_samples * sample_fraction < 1) {
    throw std::runtime_error("sample_fraction too small, no observations sampled.");
  }

  size_t num_types = observables.size();
  std::vector<std::vector<double>> observations_by_type(num_types);
  for (auto it : observables) {
    size_t type = it.first;
    size_t index = it.second;

    observations_by_type[type].resize(num_samples);
    for (size_t row = 0; row < data->get_num_rows(); ++row) {
      observations_by_type[type][row] = data->get(row, index);
    }
  }
  Observations observations(observations_by_type, data->get_num_rows());

  uint num_groups = (uint) num_trees / ci_group_size;

  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, num_groups - 1, num_threads);

  std::vector<std::future<std::vector<std::shared_ptr<Tree>>>> futures;
  futures.reserve(thread_ranges.size());

  std::vector<std::shared_ptr<Tree>> trees;
  trees.reserve(num_trees);

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_trees_batch = thread_ranges[i + 1] - start_index;
    futures.push_back(std::async(std::launch::async,
                                 &ForestTrainer::train_batch,
                                 this,
                                 start_index,
                                 num_trees_batch,
                                 data,
                                 observations));
  }

  for (auto& future : futures) {
    std::vector<std::shared_ptr<Tree>> thread_trees = future.get();
    trees.insert(trees.end(), thread_trees.begin(), thread_trees.end());
  }

  return Forest::create(trees, data, observables);
}

std::vector<std::shared_ptr<Tree>> ForestTrainer::train_batch(
    size_t start,
    size_t num_trees,
    Data* data,
    const Observations& observations) {
  std::mt19937_64 random_number_generator(seed + start);
  std::uniform_int_distribution<uint> udist;
  std::vector<std::shared_ptr<Tree>> trees;

  if (ci_group_size == 1) {
    trees.reserve(num_trees);
  } else {
    trees.reserve(num_trees * ci_group_size);
  }

  for (size_t i = 0; i < num_trees; i++) {
    uint tree_seed = udist(random_number_generator);
    SamplingOptions sampling_options(sample_with_replacement, sample_weights);
    RandomSampler sampler(tree_seed, sampling_options);

    if (ci_group_size == 1) {
      std::vector<size_t> samples;
      std::vector<size_t> oob_samples;
      sampler.sample(data->get_num_rows(), sample_fraction, samples, oob_samples);

      std::shared_ptr<Tree> tree = tree_trainer->train(data, observations,
          sampler, samples);
      tree->set_oob_samples(oob_samples);
      trees.push_back(tree);
    } else {
      std::vector<std::shared_ptr<Tree>> group = train_ci_group(
          data, observations, sampler, sample_fraction);
      trees.insert(trees.end(), group.begin(), group.end());
    }
  }
  return trees;
}

std::vector<std::shared_ptr<Tree>> ForestTrainer::train_ci_group(Data* data,
                                                                 const Observations& observations,
                                                                 RandomSampler& sampler,
                                                                 double sample_fraction) {
  std::vector<std::shared_ptr<Tree>> trees;

  std::vector<size_t> sample;
  std::vector<size_t> oob_sample;
  sampler.sample(data->get_num_rows(), 0.5, sample, oob_sample);

  for (size_t i = 0; i < ci_group_size; ++i) {
    std::vector<size_t> subsample;
    std::vector<size_t> oob_subsample;
    sampler.subsample(sample, sample_fraction * 2, subsample, oob_subsample);
    oob_subsample.insert(oob_subsample.end(), oob_sample.begin(), oob_sample.end());

    std::shared_ptr<Tree> tree = tree_trainer->train(data, observations,
        sampler, subsample);
    tree->set_oob_samples(oob_subsample);
    trees.push_back(tree);
  }

  return trees;
}

void ForestTrainer::set_split_select_weights(std::vector<double>& split_select_weights,
                                             size_t num_independent_variables) {
// Reserve space
  this->split_select_weights.resize(num_independent_variables);
  this->split_select_vars.resize(num_independent_variables);
  deterministic_vars.reserve(num_independent_variables);

  // Split up in deterministic and weighted variables, ignore zero weights
  // Size should be 1 x num_independent_variables or num_trees x num_independent_variables
  if (split_select_weights.size() != num_independent_variables) {
    throw std::runtime_error("Number of split select weights not equal to number of independent variables.");
  }

  for (size_t j = 0; j < split_select_weights.size(); ++j) {
    double weight = split_select_weights[j];

    size_t var = j;
    for (auto& skip : no_split_variables) {
      if (var >= skip) {
        ++var;
      }
    }

    if (weight == 1) {
      deterministic_vars.push_back(var);
    } else if (weight < 1 && weight > 0) {
      this->split_select_vars[j] = var;
      this->split_select_weights[j] = weight;
    } else if (weight < 0 || weight > 1) {
      throw std::runtime_error("One or more split select weights not in range [0,1].");
    }
  }
}
