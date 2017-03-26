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
                             std::shared_ptr<PredictionStrategy> prediction_strategy) :
    verbose_out(0), num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), seed(0),
    prediction_mode(false), sample_with_replacement(
    true), memory_saving_splitting(false), sample_fraction(
    1), num_threads(DEFAULT_NUM_THREADS), observables(observables),
    relabeling_strategy(relabeling_strategy),
    splitting_rule_factory(splitting_rule_factory),
    prediction_strategy(prediction_strategy) {}

void ForestTrainer::init(uint mtry,
                         uint num_trees,
                         std::ostream *verbose_out,
                         uint seed,
                         uint num_threads,
                         std::string load_forest_filename,
                         uint min_node_size,
                         std::vector<size_t> no_split_variables,
                         std::string split_select_weights_file,
                         std::vector<std::string>& always_split_variable_names,
                         bool sample_with_replacement,
                         bool memory_saving_splitting,
                         std::string case_weights_file,
                         double sample_fraction,
                         bool honesty,
                         uint ci_group_size) {

  this->verbose_out = verbose_out;
  this->always_split_variable_names = always_split_variable_names;
  this->split_select_weights_file = split_select_weights_file;
  this->case_weights_file = case_weights_file;

  // Set prediction mode
  bool prediction_mode = false;
  if (!load_forest_filename.empty()) {
    prediction_mode = true;
  }

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
  this->prediction_mode = prediction_mode;
  this->sample_with_replacement = sample_with_replacement;
  this->memory_saving_splitting = memory_saving_splitting;

  if (ci_group_size > 1 && sample_fraction > 0.5) {
    throw std::runtime_error("When confidence intervals are enabled, the"
        " sampling fraction must be less than 0.5.");
  }
  this->sample_fraction = sample_fraction;

  this->no_split_variables = no_split_variables;
  for (auto it : observables) {
    this->no_split_variables.push_back(it.second);
  }

  if (seed != 0) {
    this->seed = seed;
  } else {
    std::random_device random_device;
    this->seed = random_device();
  }

  this->ci_group_size = ci_group_size;

  // Sort no split variables in ascending order
  std::sort(no_split_variables.begin(), no_split_variables.end());

  TreeOptions tree_options(mtry,
      min_node_size,
      split_select_weights,
      split_select_varIDs,
      deterministic_varIDs,
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
  if (!case_weights_file.empty()) {
    read_vector_from_file(case_weights, case_weights_file);
    if (case_weights.size() != num_samples - 1) {
      throw std::runtime_error("Number of case weights is not equal to number of samples.");
    }
  }

  if (mtry == 0) {
    mtry = ceil((num_variables - no_split_variables.size()) / 3.0);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }

  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    set_always_split_variables(data, always_split_variable_names, num_independent_variables);
  }

  // Check if any observations samples
  if ((size_t) num_samples * sample_fraction < 1) {
    throw std::runtime_error("sample_fraction too small, no observations sampled.");
  }

  // Check if mtry is in valid range
  if (this->mtry > num_variables - 1) {
    throw std::runtime_error("mtry can not be larger than number of variables in data.");
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
  split_sequence(thread_ranges, 0, num_groups - 1, num_threads);

  std::vector<std::future<std::vector<std::shared_ptr<Tree>>>> futures;
  futures.reserve(num_threads);

  std::vector<std::shared_ptr<Tree>> trees;
  trees.reserve(num_trees);

  for (uint i = 0; i < num_threads; ++i) {
    futures.push_back(std::async(std::launch::async,
                                 &ForestTrainer::train_batch,
                                 this,
                                 i,
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
    uint thread_idx,
    Data* data,
    const Observations& observations) {
  std::mt19937_64 random_number_generator(seed + thread_idx);
  std::uniform_int_distribution<uint> udist;
  std::vector<std::shared_ptr<Tree>> trees;

  if (thread_ranges.size() > thread_idx + 1) {
    size_t start = thread_ranges[thread_idx];
    size_t end = thread_ranges[thread_idx + 1];
    if (ci_group_size == 1) {
      trees.reserve(end - start);
    } else {
      trees.reserve((end - start) * ci_group_size);
    }

    for (size_t i = start; i < end; ++i) {
      uint tree_seed = udist(random_number_generator);
      SamplingOptions sampling_options(sample_with_replacement, case_weights);
      BootstrapSampler bootstrap_sampler(tree_seed, sampling_options);

      if (ci_group_size == 1) {
        std::vector<size_t> sampleIDs;
        std::vector<size_t> oob_sampleIDs;
        bootstrap_sampler.sample(data->get_num_rows(), sample_fraction, sampleIDs, oob_sampleIDs);

        std::shared_ptr<Tree> tree = tree_trainer->train(data, observations,
            bootstrap_sampler, sampleIDs);
        tree->set_oob_sampleIDs(oob_sampleIDs);
        trees.push_back(tree);
      } else {
        std::vector<std::shared_ptr<Tree>> group = train_ci_group(
            data, observations, bootstrap_sampler, sample_fraction);
        trees.insert(trees.end(), group.begin(), group.end());
      }
    }
  }
  return trees;
}

std::vector<std::shared_ptr<Tree>> ForestTrainer::train_ci_group(Data* data,
                                                                 const Observations& observations,
                                                                 BootstrapSampler& bootstrap_sampler,
                                                                 double sample_fraction) {
  std::vector<std::shared_ptr<Tree>> trees;

  std::vector<size_t> sampleIDs;
  std::vector<size_t> oob_sampleIDs;
  bootstrap_sampler.sample(data->get_num_rows(), 0.5, sampleIDs, oob_sampleIDs);

  for (size_t i = 0; i < ci_group_size; ++i) {
    std::vector<size_t> subsampleIDs;
    std::vector<size_t> oob_subsampleIDs;
    bootstrap_sampler.subsample(sampleIDs, sample_fraction * 2, subsampleIDs, oob_subsampleIDs);
    oob_subsampleIDs.insert(oob_subsampleIDs.end(), oob_sampleIDs.begin(), oob_sampleIDs.end());

    std::shared_ptr<Tree> tree = tree_trainer->train(data, observations,
        bootstrap_sampler, subsampleIDs);
    tree->set_oob_sampleIDs(oob_subsampleIDs);
    trees.push_back(tree);
  }

  return trees;
}

void ForestTrainer::set_split_select_weights(std::vector<double>& split_select_weights,
                                             size_t num_independent_variables) {
// Reserve space
  this->split_select_weights.resize(num_independent_variables);
  this->split_select_varIDs.resize(num_independent_variables);
  deterministic_varIDs.reserve(num_independent_variables);

  // Split up in deterministic and weighted variables, ignore zero weights
  // Size should be 1 x num_independent_variables or num_trees x num_independent_variables
  if (split_select_weights.size() != num_independent_variables) {
    throw std::runtime_error("Number of split select weights not equal to number of independent variables.");
  }

  for (size_t j = 0; j < split_select_weights.size(); ++j) {
    double weight = split_select_weights[j];

    size_t varID = j;
    for (auto& skip : no_split_variables) {
      if (varID >= skip) {
        ++varID;
      }
    }

    if (weight == 1) {
      deterministic_varIDs.push_back(varID);
    } else if (weight < 1 && weight > 0) {
      this->split_select_varIDs[j] = varID;
      this->split_select_weights[j] = weight;
    } else if (weight < 0 || weight > 1) {
      throw std::runtime_error("One or more split select weights not in range [0,1].");
    }
  }

  if (deterministic_varIDs.size() > this->mtry) {
    throw std::runtime_error("Number of ones in split select weights cannot be larger than mtry.");
  }
  if (deterministic_varIDs.size() + split_select_varIDs.size() < mtry) {
    throw std::runtime_error("Too many zeros in split select weights. Need at least mtry variables to split at.");
  }
}

void ForestTrainer::set_always_split_variables(Data* data,
                                               std::vector<std::string>& always_split_variable_names,
                                               size_t num_independent_variables) {
  deterministic_varIDs.clear();
  deterministic_varIDs.reserve(num_independent_variables);

  for (auto& variable_name : always_split_variable_names) {
    size_t varID = data->get_variable_id(variable_name);
    deterministic_varIDs.push_back(varID);
  }

  if (deterministic_varIDs.size() + this->mtry > num_independent_variables) {
    throw std::runtime_error(
        "Number of variables to be always considered for splitting plus mtry cannot be larger than number of independent variables.");
  }
}
