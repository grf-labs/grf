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
#include <math.h>
#ifndef WIN_R_BUILD
#include <thread>
#include <chrono>
#endif

#include "utility.h"
#include "Forest.h"
#include "DataChar.h"
#include "DataDouble.h"
#include "DataFloat.h"

Forest::Forest() :
    verbose_out(0), num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), num_variables(0), num_independent_variables(
        0), seed(0), dependent_varID(0), num_samples(0), prediction_mode(false), memory_mode(MEM_DOUBLE), sample_with_replacement(
        true), num_threads(DEFAULT_NUM_THREADS), data(0), overall_prediction_error(0), importance_mode(
        DEFAULT_IMPORTANCE_MODE), progress(0) {
}

Forest::~Forest() {
  for (auto& tree : trees) {
    delete tree;
  }
}

void Forest::initCpp(std::string dependent_variable_name, MemoryMode memory_mode, std::string input_file, uint mtry,
    std::string output_prefix, uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
    std::string load_forest_filename, ImportanceMode importance_mode, uint min_node_size,
    std::string split_select_weights_file, std::vector<std::string>& always_split_variable_names,
    std::string status_variable_name, bool sample_with_replacement,
    std::vector<std::string>& unordered_variable_names) {

  this->verbose_out = verbose_out;

  // Initialize data with memmode
  switch (memory_mode) {
  case MEM_DOUBLE:
    data = new DataDouble();
    break;
  case MEM_FLOAT:
    data = new DataFloat();
    break;
  case MEM_CHAR:
    data = new DataChar();
    break;
  }

  // Load data
  *verbose_out << "Loading input file: " << input_file << "." << std::endl;
  bool rounding_error = data->loadFromFile(input_file);
  if (rounding_error) {
    *verbose_out << "Warning: Rounding or Integer overflow occurred. Use FLOAT or DOUBLE precision to avoid this."
        << std::endl;
  }

  // Set prediction mode
  bool prediction_mode = false;
  if (!load_forest_filename.empty()) {
    prediction_mode = true;
  }

  // Call other init function
  init(dependent_variable_name, memory_mode, data, mtry, output_prefix, num_trees, seed, num_threads, importance_mode,
      min_node_size, status_variable_name, prediction_mode, sample_with_replacement, unordered_variable_names);

  if (prediction_mode) {
    loadFromFile(load_forest_filename);
  }
  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    setAlwaysSplitVariables(always_split_variable_names);
  }

  // Load split select weights from file
  if (!split_select_weights_file.empty()) {
    std::vector<double> split_select_weights;
    loadDoubleVectorFromFile(split_select_weights, split_select_weights_file);
    if (split_select_weights.size() != num_variables - 1) {
      throw std::runtime_error("Number of split select weights is not equal to number of independent variables.");
    }
    setSplitWeightVector(split_select_weights);
  }

  // Check if all catvars are coded in integers starting at 1
  if (!unordered_variable_names.empty()) {
    std::string error_message = checkUnorderedVariables(data, unordered_variable_names);
    if (!error_message.empty()) {
      throw std::runtime_error(error_message);
    }
  }
}

void Forest::initR(std::string dependent_variable_name, Data* input_data, uint mtry, uint num_trees,
    std::ostream* verbose_out, uint seed, uint num_threads, ImportanceMode importance_mode, uint min_node_size,
    std::vector<double>& split_select_weights, std::vector<std::string>& always_split_variable_names,
    std::string status_variable_name, bool prediction_mode, bool sample_with_replacement,
    std::vector<std::string>& unordered_variable_names) {

  this->verbose_out = verbose_out;

  // Call other init function
  init(dependent_variable_name, MEM_DOUBLE, input_data, mtry, "", num_trees, seed, num_threads, importance_mode,
      min_node_size, status_variable_name, prediction_mode, sample_with_replacement, unordered_variable_names);

  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    setAlwaysSplitVariables(always_split_variable_names);
  }

  // Set split select weights
  if (!split_select_weights.empty()) {
    setSplitWeightVector(split_select_weights);
  }
}

void Forest::init(std::string dependent_variable_name, MemoryMode memory_mode, Data* input_data, uint mtry,
    std::string output_prefix, uint num_trees, uint seed, uint num_threads, ImportanceMode importance_mode,
    uint min_node_size, std::string status_variable_name, bool prediction_mode, bool sample_with_replacement,
    std::vector<std::string>& unordered_variable_names) {

  // Initialize data with memmode
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
#ifdef WIN_R_BUILD
    this->num_threads = 1;
#else
    this->num_threads = std::thread::hardware_concurrency();
#endif
  } else {
    this->num_threads = num_threads;
  }

  // Set member variables
  this->num_trees = num_trees;
  this->mtry = mtry;
  this->seed = seed;
  this->output_prefix = output_prefix;
  this->importance_mode = importance_mode;
  this->min_node_size = min_node_size;
  this->memory_mode = memory_mode;
  this->prediction_mode = prediction_mode;
  this->sample_with_replacement = sample_with_replacement;

  // Set number of samples and variables
  num_samples = data->getNumRows();
  num_variables = data->getNumCols();

  // Convert dependent variable name to ID
  if (!prediction_mode && !dependent_variable_name.empty()) {
    dependent_varID = data->getVariableID(dependent_variable_name);
  }

  // Set unordered factor variables
  if (!prediction_mode) {
    is_ordered_variable.resize(num_variables, true);
    for (auto& variable_name : unordered_variable_names) {
      size_t varID = data->getVariableID(variable_name);
      is_ordered_variable[varID] = false;
    }
  }

  no_split_variables.push_back(dependent_varID);

  initInternal(status_variable_name);

  num_independent_variables = num_variables - no_split_variables.size();

  // Sort no split variables in ascending order
  std::sort(no_split_variables.begin(), no_split_variables.end());

  // Check if mtry is in valid range
  if (this->mtry > num_variables - 1) {
    throw std::runtime_error("mtry can not be larger than number of variables in data.");
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
    }
    computePredictionError();

    if (importance_mode == IMP_GINI) {
      if (verbose) {
        *verbose_out << "Computing variable importance .." << std::endl;
      }
      computeGiniImportance();
    } else if (importance_mode > IMP_GINI) {
      if (verbose) {
        *verbose_out << "Computing permutation variable importance .." << std::endl;
      }
      computePermutationImportance();
    }
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
  *verbose_out << "Variable importance mode:          " << importance_mode << std::endl;
  *verbose_out << "Memory mode:                       " << memory_mode << std::endl;
  *verbose_out << "Seed:                              " << seed << std::endl;
  *verbose_out << "Number of threads:                 " << num_threads << std::endl;
  *verbose_out << std::endl;

  if (prediction_mode) {
    writePredictionFile();
  } else {
    *verbose_out << "Overall OOB prediction error:      " << overall_prediction_error << std::endl;
    *verbose_out << std::endl;

    if (!split_select_weights.empty()) {
      *verbose_out
          << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
          << std::endl;
    }

    if (importance_mode != IMP_NONE) {
      writeImportanceFile();
    }

    writeConfusionFile();
  }
}

void Forest::writeImportanceFile() {

  // Open importance file for writing
  std::string filename = output_prefix + ".importance";
  std::ofstream importance_file;
  importance_file.open(filename, std::ios::out);
  if (!importance_file.good()) {
    throw std::runtime_error("Could not write to importance file: " + filename + ".");
  }

  // Write importance to file
  for (size_t i = 0; i < variable_importance.size(); ++i) {
    size_t varID = i;
    for (auto& skip : no_split_variables) {
      if (varID >= skip) {
        ++varID;
      }
    }
    std::string variable_name = data->getVariableNames()[varID];
    importance_file << variable_name << ": " << variable_importance[i] << std::endl;
  }

  importance_file.close();
  *verbose_out << "Saved variable importance to file " << filename << "." << std::endl;
}

void Forest::saveToFile() {

  // Open file for writing
  std::string filename = output_prefix + ".forest";
  std::ofstream outfile;
  outfile.open(filename, std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: " + filename + ".");
  }

  // Write dependent_varID
  outfile.write((char*) &dependent_varID, sizeof(dependent_varID));

  // Write num_trees
  outfile.write((char*) &num_trees, sizeof(num_trees));

  // Write is_ordered_variable
  saveVector1D(is_ordered_variable, outfile);

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
    trees[i]->init(data, mtry, dependent_varID, num_samples, tree_seed, &deterministic_varIDs, &split_select_varIDs,
        &split_select_weights, importance_mode, min_node_size, &no_split_variables, sample_with_replacement,
        &is_ordered_variable);
  }

  // Grow trees in multiple threads
#ifdef WIN_R_BUILD
  progress = 0;
  clock_t start_time = clock();
  clock_t lap_time = clock();
  for (size_t i = 0; i < num_trees; ++i) {
    trees[i]->grow();
    progress++;
    showProgress("Growing trees..", start_time, lap_time);
  }
#else
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
#endif

}

void Forest::predict() {

  // Predict trees in multiple threads and join the threads with the main thread
#ifdef WIN_R_BUILD
  progress = 0;
  clock_t start_time = clock();
  clock_t lap_time = clock();
  for (size_t i = 0; i < num_trees; ++i) {
    trees[i]->predict(data, false);
    progress++;
    showProgress("Predicting..", start_time, lap_time);
  }
#else
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
#endif

  // Call special functions for subclasses
  predictInternal();
}

void Forest::computePredictionError() {

  // Predict trees in multiple threads
#ifdef WIN_R_BUILD
  progress = 0;
  clock_t start_time = clock();
  clock_t lap_time = clock();
  for (size_t i = 0; i < num_trees; ++i) {
    trees[i]->predict(data, true);
    progress++;
    showProgress("Predicting..", start_time, lap_time);
  }
#else
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (uint i = 0; i < num_threads; ++i) {
    threads.push_back(std::thread(&Forest::predictTreesInThread, this, i, data, true));
  }
  for (auto &thread : threads) {
    thread.join();
  }
#endif

  // Call special function for subclasses
  computePredictionErrorInternal();
}

void Forest::computeGiniImportance() {

  // Initialize with 0.
  variable_importance.resize(num_independent_variables, 0);

  // Sum tree importance and divide by number of trees
  for (size_t t = 0; t < trees.size(); ++t) {
    std::vector<double> tree_importance = trees[t]->getVariableImportance();
    for (size_t i = 0; i < tree_importance.size(); ++i) {
      variable_importance[i] += tree_importance[i];
    }
  }
  for (auto& v : variable_importance) {
    v /= num_trees;
  }
}

void Forest::computePermutationImportance() {

  // Compute tree permutation importance in multiple threads
#ifdef WIN_R_BUILD
  progress = 0;
  clock_t start_time = clock();
  clock_t lap_time = clock();
  for (size_t i = 0; i < num_trees; ++i) {
    trees[i]->computePermutationImportance();
    progress++;
    showProgress("Computing permutation importance..", start_time, lap_time);
  }
#else
  progress = 0;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (uint i = 0; i < num_threads; ++i) {
    threads.push_back(std::thread(&Forest::computeTreePermutationImportanceInThread, this, i));
  }
  showProgress("Computing permutation importance..");
  for (auto &thread : threads) {
    thread.join();
  }
#endif

  // Initialize with 0.
  std::vector<double> variance;
  variance.resize(num_independent_variables, 0);
  variable_importance.resize(num_independent_variables, 0);

  // Sum tree importance and divide by number of trees
  for (size_t t = 0; t < trees.size(); ++t) {
    std::vector<double> tree_importance = trees[t]->getVariableImportance();
    for (size_t i = 0; i < tree_importance.size(); ++i) {
      variable_importance[i] += tree_importance[i];

      // Variance for scaled permutation importance
      if (importance_mode == IMP_PERM_BREIMAN) {
        variance[i] += tree_importance[i] * tree_importance[i];
      } else if (importance_mode == IMP_PERM_LIAW) {
        variance[i] += tree_importance[i] * tree_importance[i] * trees[i]->getNumSamplesOob();
      }
    }
  }
  for (size_t i = 0; i < variable_importance.size(); ++i) {
    variable_importance[i] /= num_trees;

    // Normalize by variance for scaled permutation importance
    if (importance_mode == IMP_PERM_BREIMAN || importance_mode == IMP_PERM_LIAW) {
      if (variance[i] != 0) {
        variance[i] = variance[i] / num_trees - variable_importance[i] * variable_importance[i];
        variable_importance[i] /= sqrt(variance[i] / num_trees);
      }
    }
  }
}

#ifndef WIN_R_BUILD
void Forest::growTreesInThread(uint thread_idx) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->grow();

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
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

void Forest::computeTreePermutationImportanceInThread(uint thread_idx) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      trees[i]->computePermutationImportance();

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}
#endif

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

  // Read is_ordered_variable
  readVector1D(is_ordered_variable, infile);

  // Read tree data. This is different for tree types -> virtual function
  loadFromFileInternal(infile);

  infile.close();

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void Forest::setSplitWeightVector(std::vector<double>& split_select_weights) {

  if (split_select_weights.size() != num_independent_variables) {
    throw std::runtime_error("Number of split select weights not equal to number of independent variables.");
  }

  // Split up in deterministic and weighted variables, ignore zero weights
  for (size_t i = 0; i < split_select_weights.size(); ++i) {
    double weight = split_select_weights[i];
    size_t varID = i;
    for (auto& skip : no_split_variables) {
      if (varID >= skip) {
        ++varID;
      }
    }

    if (weight == 1) {
      deterministic_varIDs.push_back(varID);
    } else if (weight < 1 && weight > 0) {
      this->split_select_varIDs.push_back(varID);
      this->split_select_weights.push_back(weight);
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

void Forest::setAlwaysSplitVariables(std::vector<std::string>& always_split_variable_names) {

  for (auto& variable_name : always_split_variable_names) {
    size_t varID = data->getVariableID(variable_name);
    deterministic_varIDs.push_back(varID);
  }

  if (deterministic_varIDs.size() + this->mtry > num_independent_variables) {
    throw std::runtime_error(
        "Number of variables to be always considered for splitting plus mtry cannot be larger than number of independent variables.");
  }
}

#ifdef WIN_R_BUILD
void Forest::showProgress(std::string operation, clock_t start_time, clock_t& lap_time) {

  double elapsed_time = (clock() - lap_time) / CLOCKS_PER_SEC;
  if (elapsed_time > STATUS_INTERVAL) {
    double relative_progress = (double) progress / (double) num_trees;
    double time_from_start = (clock() - start_time) / CLOCKS_PER_SEC;
    uint remaining_time = (1 / relative_progress - 1) * time_from_start;
    *verbose_out << operation << " Progress: " << round(100 * relative_progress) << "%. Estimated remaining time: "
    << beautifyTime(remaining_time) << "." << std::endl;
    lap_time = clock();
  }
}
#else
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
#endif
