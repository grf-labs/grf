#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <thread>
#include <future>
#include "utility.h"
#include "ForestTrainer.h"

ForestTrainer::ForestTrainer(std::unordered_map<std::string, size_t> observables,
                         RelabelingStrategy *relabeling_strategy,
                         SplittingRule *splitting_rule,
                         PredictionStrategy *prediction_strategy) :
    verbose_out(0), num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), seed(0),
    prediction_mode(false), sample_with_replacement(
    true), memory_saving_splitting(false), keep_inbag(false), sample_fraction(
    1), num_threads(DEFAULT_NUM_THREADS), observables(observables),
    relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule), prediction_strategy(prediction_strategy) {
}

void ForestTrainer::init(uint mtry,
                       uint num_trees, std::ostream *verbose_out, uint seed, uint num_threads,
                       std::string load_forest_filename, uint min_node_size,
                       std::string split_select_weights_file, std::vector<std::string> &always_split_variable_names,
                       bool sample_with_replacement, bool memory_saving_splitting,
                       std::string case_weights_file,
                       double sample_fraction) {

  this->verbose_out = verbose_out;
  this->always_split_variable_names = always_split_variable_names;
  this->split_select_weights_file = split_select_weights_file;
  this->case_weights_file = case_weights_file;

  // Set prediction mode
  bool prediction_mode = false;
  if (!load_forest_filename.empty()) {
    prediction_mode = true;
  }

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
  this->prediction_mode = prediction_mode;
  this->sample_with_replacement = sample_with_replacement;
  this->memory_saving_splitting = memory_saving_splitting;
  this->sample_fraction = sample_fraction;


  for (auto it : observables) {
    no_split_variables.push_back(it.second);
  }

  // Sort no split variables in ascending order
  std::sort(no_split_variables.begin(), no_split_variables.end());

  TreeOptions* tree_options = new TreeOptions(
      mtry,
      min_node_size,
      &split_select_weights,
      &split_select_varIDs,
      &deterministic_varIDs,
      &no_split_variables);
  tree_model = new TreeModel(relabeling_strategy,
                             splitting_rule,
                             prediction_strategy,
                             tree_options);

}

Forest* ForestTrainer::train(Data* data) {
  size_t num_samples = data->getNumRows();
  size_t num_variables = data->getNumCols();
  size_t num_independent_variables = num_variables - no_split_variables.size();

  if (!split_select_weights_file.empty()) {
    std::vector<double> split_select_weights;
    split_select_weights.resize(1);
    loadDoubleVectorFromFile(split_select_weights, split_select_weights_file);
    if (split_select_weights.size() != num_variables - 1) {
      throw std::runtime_error("Number of split select weights is not equal to number of independent variables.");
    }
    setSplitWeightVector(split_select_weights, num_independent_variables);
  }

  // Load case weights from file
  if (!case_weights_file.empty()) {
    loadDoubleVectorFromFile(case_weights, case_weights_file);
    if (case_weights.size() != num_samples - 1) {
      throw std::runtime_error("Number of case weights is not equal to number of samples.");
    }
  }

  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - no_split_variables.size()));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }

  // Set variables to be always considered for splitting
  if (!always_split_variable_names.empty()) {
    setAlwaysSplitVariables(data, always_split_variable_names, num_independent_variables);
  }

  // Check if any observations samples
  if ((size_t) num_samples * sample_fraction < 1) {
    throw std::runtime_error("sample_fraction too small, no observations sampled.");
  }

  // Check if mtry is in valid range
  if (this->mtry > num_variables - 1) {
    throw std::runtime_error("mtry can not be larger than number of variables in data.");
  }

  auto observations_by_type = new std::unordered_map<std::string, std::vector<double>>();
  for (auto it : observables) {
    std::string name = it.first;
    size_t index = it.second;

    for (int row = 0; row < data->getNumRows(); row++) {
      (*observations_by_type)[name].push_back(data->get(row, index));
    }
  }
  Observations* observations = new Observations(*observations_by_type, data->getNumRows());

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::vector<std::future<std::vector<Tree*>>> futures;
  futures.reserve(num_threads);

  std::vector<Tree*>* trees = new std::vector<Tree*>();
  trees->reserve(num_trees);

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::vector<Tree*>> promise;
    std::future<std::vector<Tree*>> future = promise.get_future();
    threads.push_back(std::thread(&ForestTrainer::growTreesInThread,
                                  this,
                                  i,
                                  data,
                                  observations,
                                  std::move(promise)));
    futures.push_back(std::move(future));
  }

  for (auto &future : futures) {
    future.wait();
    std::vector<Tree*> thread_trees = future.get();
    trees->insert(trees->end(), thread_trees.begin(), thread_trees.end());
  }

  for (auto& thread : threads) {
    thread.join();
  }

  return new Forest(trees, data, observables);
}

void ForestTrainer::growTreesInThread(uint thread_idx,
                                    Data* data,
                                    Observations* observations,
                                    std::promise<std::vector<Tree *>> promise) {
  std::uniform_int_distribution<uint> udist;
  std::vector<Tree*> trees;
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      uint tree_seed;
      if (seed == 0) {
        tree_seed = udist(random_number_generator);
      } else {
        tree_seed = (i + 1) * seed;
      }

      SamplingOptions* sampling_options = new SamplingOptions(sample_with_replacement,
                                                              sample_fraction,
                                                              &case_weights,
                                                              keep_inbag);
      BootstrapSampler* bootstrap_sampler = new BootstrapSampler(data->getNumRows(),
                                                                 tree_seed,
                                                                 sampling_options);

      trees.push_back(tree_model->train(data,
                                        bootstrap_sampler,
                                        observations));
    }
  }
  promise.set_value(trees);
}

void ForestTrainer::setSplitWeightVector(std::vector<double> &split_select_weights,
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
    for (auto &skip : no_split_variables) {
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

void ForestTrainer::setAlwaysSplitVariables(Data* data,
                                          std::vector<std::string>& always_split_variable_names,
                                          size_t num_independent_variables) {
  deterministic_varIDs.clear();
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

void ForestTrainer::writeConfusionFile(Data* prediction_data, std::vector<std::vector<double>> predictions) {

// Open confusion file for writing
  std::string filename = "gradientforest.confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  // Write
  outfile << "Prediction X1" << std::endl;
  for (size_t i = 0; i < predictions.size(); ++i) {
    for (size_t j = 0; j < predictions[i].size(); ++j) {
      outfile << predictions[i][j] << " " << prediction_data->get(i, 0) << " ";
    }
    outfile << std::endl;
  }

  outfile.close();
  *verbose_out << "Saved prediction error to file " << filename << "." << std::endl;
}
