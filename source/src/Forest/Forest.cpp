/*-------------------------------------------------------------------------------
 This file is part of Ranger.

 Ranger is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Ranger is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Ranger. If not, see <http://www.gnu.org/licenses/>.

 Written by:

 Marvin N. Wright
 Institut für Medizinische Biometrie und Statistik
 Universität zu Lübeck
 Ratzeburger Allee 160
 23562 Lübeck

 http://www.imbs-luebeck.de
 wright@imbs.uni-luebeck.de
 #-------------------------------------------------------------------------------*/

#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <thread>
#include <chrono>
#include "utility.h"
#include "Forest.h"
#include "DataChar.h"
#include "DataDouble.h"
#include "DataFloat.h"

Forest::Forest(RelabelingStrategy* relabeling_strategy,
               SplittingRule* splitting_rule) :
    verbose_out(0), num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), num_variables(0), num_independent_variables(
        0), seed(0), dependent_varID(0), num_samples(0), prediction_mode(false), memory_mode(MEM_DOUBLE), sample_with_replacement(
        true), memory_saving_splitting(false), keep_inbag(false), sample_fraction(
        1), num_threads(DEFAULT_NUM_THREADS), data(0), overall_prediction_error(
        0), progress(0), relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule) {
}

Forest::~Forest() {
  for (auto& tree : trees) {
    delete tree;
  }
}

void Forest::initCpp(std::string dependent_variable_name, MemoryMode memory_mode, uint mtry,
    uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
    std::string load_forest_filename, uint min_node_size,
    std::string split_select_weights_file, std::vector<std::string>& always_split_variable_names,
    std::string status_variable_name, bool sample_with_replacement, bool memory_saving_splitting,
    std::string case_weights_file,
    double sample_fraction, Data* data) {

  this->verbose_out = verbose_out;

  // Set prediction mode
  bool prediction_mode = false;
  if (!load_forest_filename.empty()) {
    prediction_mode = true;
  }

  // Call other init function
  init(dependent_variable_name, memory_mode, data, mtry, num_trees, seed, num_threads,
      min_node_size, status_variable_name, prediction_mode, sample_with_replacement,
      memory_saving_splitting, sample_fraction);

  if (prediction_mode) {
    loadFromFile(load_forest_filename);
  }
  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    setAlwaysSplitVariables(always_split_variable_names);
  }

  // TODO: Read 2d weights for tree-wise split select weights
  // Load split select weights from file
  if (!split_select_weights_file.empty()) {
    std::vector<std::vector<double>> split_select_weights;
    split_select_weights.resize(1);
    loadDoubleVectorFromFile(split_select_weights[0], split_select_weights_file);
    if (split_select_weights[0].size() != num_variables - 1) {
      throw std::runtime_error("Number of split select weights is not equal to number of independent variables.");
    }
    setSplitWeightVector(split_select_weights);
  }

  // Load case weights from file
  if (!case_weights_file.empty()) {
    loadDoubleVectorFromFile(case_weights, case_weights_file);
    if (case_weights.size() != num_samples - 1) {
      throw std::runtime_error("Number of case weights is not equal to number of samples.");
    }
  }
}

void Forest::init(std::string dependent_variable_name, MemoryMode memory_mode, Data* input_data, uint mtry,
    uint num_trees, uint seed, uint num_threads,
    uint min_node_size, std::string status_variable_name, bool prediction_mode, bool sample_with_replacement,
    bool memory_saving_splitting,
    double sample_fraction) {

  this->data = input_data;

  // Initialize random number generator and set seed
  if (seed == 0) {
    std::random_device random_device;
    random_number_generator.seed(random_device());
  } else {
    random_number_generator.seed(seed);
  }

  // Set number of threads
  if (num_threads == DEFAULT_NUM_THREADS) {
    this->num_threads = std::thread::hardware_concurrency();
  } else {
    this->num_threads = num_threads;
  }

  // Set member variables
  this->num_trees = num_trees;
  this->mtry = mtry;
  this->seed = seed;
  this->min_node_size = min_node_size;
  this->memory_mode = memory_mode;
  this->prediction_mode = prediction_mode;
  this->sample_with_replacement = sample_with_replacement;
  this->memory_saving_splitting = memory_saving_splitting;
  this->sample_fraction = sample_fraction;

  // Set number of samples and variables
  num_samples = data->getNumRows();
  num_variables = data->getNumCols();

  // Convert dependent variable name to ID
  if (!prediction_mode && !dependent_variable_name.empty()) {
    dependent_varID = data->getVariableID(dependent_variable_name);
  }

  no_split_variables.push_back(dependent_varID);

  initInternal(status_variable_name);

  num_independent_variables = num_variables - no_split_variables.size();

  // Sort no split variables in ascending order
  std::sort(no_split_variables.begin(), no_split_variables.end());

  // Init split select weights
  split_select_weights.push_back(std::vector<double>());

  // Check if mtry is in valid range
  if (this->mtry > num_variables - 1) {
    throw std::runtime_error("mtry can not be larger than number of variables in data.");
  }

  // Check if any observations samples
  if ((size_t) num_samples * sample_fraction < 1) {
    throw std::runtime_error("sample_fraction too small, no observations sampled.");
  }
}

void Forest::run(bool verbose) {

  if (prediction_mode) {
    if (verbose) {
      *verbose_out << "Predicting .." << std::endl;
    }

    predict();
  } else {
    if (verbose) {
      *verbose_out << "Growing trees .." << std::endl;
    }

    grow();

    if (verbose) {
      *verbose_out << "Computing prediction error .." << std::endl;
      *verbose_out << "OOB PREDICTION OFF!!" << std::endl;
    }

    // computePredictionError();
  }
}

void Forest::writeOutput() {

  *verbose_out << std::endl;
  writeOutputInternal();
  *verbose_out << "Dependent variable name:           " << data->getVariableNames()[dependent_varID] << std::endl;
  *verbose_out << "Dependent variable ID:             " << dependent_varID << std::endl;
  *verbose_out << "Number of trees:                   " << num_trees << std::endl;
  *verbose_out << "Sample size:                       " << num_samples << std::endl;
  *verbose_out << "Number of independent variables:   " << num_independent_variables << std::endl;
  *verbose_out << "Mtry:                              " << mtry << std::endl;
  *verbose_out << "Target node size:                  " << min_node_size << std::endl;
  *verbose_out << "Memory mode:                       " << memory_mode << std::endl;
  *verbose_out << "Seed:                              " << seed << std::endl;
  *verbose_out << "Number of threads:                 " << num_threads << std::endl;
  *verbose_out << std::endl;

  if (prediction_mode) {
    writePredictionFile();
  } else {
    *verbose_out << "Overall OOB prediction error:      " << overall_prediction_error << std::endl;
    *verbose_out << std::endl;

    if (!split_select_weights.empty() & !split_select_weights[0].empty()) {
      *verbose_out
          << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
          << std::endl;
    }

    writeConfusionFile();
  }
}

void Forest::saveToFile() {

  // Open file for writing
  std::string filename = "ranger.forest";
  std::ofstream outfile;
  outfile.open(filename, std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: " + filename + ".");
  }

  // Write dependent_varID
  outfile.write((char*) &dependent_varID, sizeof(dependent_varID));

  // Write num_trees
  outfile.write((char*) &num_trees, sizeof(num_trees));

  saveToFileInternal(outfile);

  // Write tree data for each tree
  for (auto& tree : trees) {
    tree->appendToFile(outfile);
  }

  // Close file
  outfile.close();
  *verbose_out << "Saved forest to file " << filename << "." << std::endl;
}

void Forest::grow() {

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);

  // Call special grow functions of subclasses. There trees must be created.
  growInternal();

  // Init trees, create a seed for each tree, based on main seed
  std::uniform_int_distribution<uint> udist;
  for (size_t i = 0; i < num_trees; ++i) {
    uint tree_seed;
    if (seed == 0) {
      tree_seed = udist(random_number_generator);
    } else {
      tree_seed = (i + 1) * seed;
    }

    // Get split select weights for tree
    std::vector<double>* tree_split_select_weights;
    if (split_select_weights.size() > 1) {
      tree_split_select_weights = &split_select_weights[i];
    } else {
      tree_split_select_weights = &split_select_weights[0];
    }

    trees[i]->init(data, mtry, dependent_varID, num_samples, tree_seed, &deterministic_varIDs, &split_select_varIDs,
        tree_split_select_weights, min_node_size, &no_split_variables, sample_with_replacement,
        &case_weights, keep_inbag, sample_fraction);
  }

  progress = 0;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (uint i = 0; i < num_threads; ++i) {
    threads.push_back(std::thread(&Forest::growTreesInThread, this, i));
  }
  showProgress("Growing trees..");
  for (auto &thread : threads) {
    thread.join();
  }
}

void Forest::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeFactory(relabeling_strategy, splitting_rule));
  }
}

void Forest::predict() {

// Predict trees in multiple threads and join the threads with the main thread
  progress = 0;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (uint i = 0; i < num_threads; ++i) {
    threads.push_back(std::thread(&Forest::predictTreesInThread, this, i, data, false));
  }
  showProgress("Predicting..");
  for (auto &thread : threads) {
    thread.join();
  }

// Call special functions for subclasses

  predictInternal();
}

void Forest::computePredictionError() {

// Predict trees in multiple threads
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (uint i = 0; i < num_threads; ++i) {
    threads.push_back(std::thread(&Forest::predictTreesInThread, this, i, data, true));
  }
  showProgress("Computing prediction error..");
  for (auto &thread : threads) {
    thread.join();
  }
  // Call special function for subclasses
  computePredictionErrorInternal();
}

void Forest::growTreesInThread(uint thread_idx) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->grow();

      // Check for user interrupt

      // Increase progress by 1 tree
      ++progress;
      condition_variable.notify_one();
    }
  }
}

void Forest::predictTreesInThread(uint thread_idx, const Data* prediction_data, bool oob_prediction) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->predict(prediction_data, oob_prediction);

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}

void Forest::loadFromFile(std::string filename) {
  *verbose_out << "Loading forest from file " << filename << "." << std::endl;

// Open file for reading
  std::ifstream infile;
  infile.open(filename, std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: " + filename + ".");
  }

// Read dependent_varID and num_trees
  infile.read((char*) &dependent_varID, sizeof(dependent_varID));
  infile.read((char*) &num_trees, sizeof(num_trees));

// Read tree data. This is different for tree types -> virtual function
  loadFromFileInternal(infile);

  infile.close();

// Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void Forest::setSplitWeightVector(std::vector<std::vector<double>>& split_select_weights) {

// Size should be 1 x num_independent_variables or num_trees x num_independent_variables
  if (split_select_weights.size() != 1 && split_select_weights.size() != num_trees) {
    throw std::runtime_error("Size of split select weights not equal to 1 or number of trees.");
  }

// Reserve space
  if (split_select_weights.size() == 1) {
    this->split_select_weights[0].resize(num_independent_variables);
  } else {
    this->split_select_weights.clear();
    this->split_select_weights.resize(num_trees, std::vector<double>(num_independent_variables));
  }
  this->split_select_varIDs.resize(num_independent_variables);
  deterministic_varIDs.reserve(num_independent_variables);

// Split up in deterministic and weighted variables, ignore zero weights
  for (size_t i = 0; i < split_select_weights.size(); ++i) {

    // Size should be 1 x num_independent_variables or num_trees x num_independent_variables
    if (split_select_weights[i].size() != num_independent_variables) {
      throw std::runtime_error("Number of split select weights not equal to number of independent variables.");
    }

    for (size_t j = 0; j < split_select_weights[i].size(); ++j) {
      double weight = split_select_weights[i][j];

      if (i == 0) {
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
          this->split_select_weights[i][j] = weight;
        } else if (weight < 0 || weight > 1) {
          throw std::runtime_error("One or more split select weights not in range [0,1].");
        }

      } else {
        if (weight < 1 && weight > 0) {
          this->split_select_weights[i][j] = weight;
        } else if (weight < 0 || weight > 1) {
          throw std::runtime_error("One or more split select weights not in range [0,1].");
        }
      }
    }
  }

  if (deterministic_varIDs.size() > this->mtry) {
    throw std::runtime_error("Number of ones in split select weights cannot be larger than mtry.");
  }
  if (deterministic_varIDs.size() + split_select_varIDs.size() < mtry) {
    throw std::runtime_error("Too many zeros in split select weights. Need at least mtry variables to split at.");
  }
}

void Forest::setAlwaysSplitVariables(std::vector<std::string>& always_split_variable_names) {

  deterministic_varIDs.reserve(num_independent_variables);

  for (auto& variable_name : always_split_variable_names) {
    size_t varID = data->getVariableID(variable_name);
    deterministic_varIDs.push_back(varID);
  }

  if (deterministic_varIDs.size() + this->mtry > num_independent_variables) {
    throw std::runtime_error(
        "Number of variables to be always considered for splitting plus mtry cannot be larger than number of independent variables.");
  }
}

void Forest::showProgress(std::string operation) {
  using std::chrono::steady_clock;
  using std::chrono::duration_cast;
  using std::chrono::seconds;

  steady_clock::time_point start_time = steady_clock::now();
  steady_clock::time_point last_time = steady_clock::now();
  std::unique_lock<std::mutex> lock(mutex);

// Wait for message from threads and show output if enough time elapsed
  while (progress < num_trees) {
    condition_variable.wait(lock);
    seconds elapsed_time = duration_cast<seconds>(steady_clock::now() - last_time);

    // Check for user interrupt
    if (progress > 0 && elapsed_time.count() > STATUS_INTERVAL) {
      double relative_progress = (double) progress / (double) num_trees;
      seconds time_from_start = duration_cast<seconds>(steady_clock::now() - start_time);
      uint remaining_time = (1 / relative_progress - 1) * time_from_start.count();
      *verbose_out << operation << " Progress: " << round(100 * relative_progress) << "%. Estimated remaining time: "
          << beautifyTime(remaining_time) << "." << std::endl;
      last_time = steady_clock::now();
    }
  }
}
