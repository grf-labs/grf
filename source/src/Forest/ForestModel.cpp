#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <thread>
#include <future>
#include "utility.h"
#include "ForestModel.h"

ForestModel::ForestModel(std::unordered_map<std::string, size_t> observables,
                         RelabelingStrategy *relabeling_strategy,
                         SplittingRule *splitting_rule,
                         PredictionStrategy *prediction_strategy) :
    verbose_out(0), num_trees(DEFAULT_NUM_TREE), mtry(0), min_node_size(0), seed(0), dependent_varID(0),
    prediction_mode(false), sample_with_replacement(
    true), memory_saving_splitting(false), keep_inbag(false), sample_fraction(
    1), num_threads(DEFAULT_NUM_THREADS), progress(0), observables(observables),
    relabeling_strategy(relabeling_strategy), splitting_rule(splitting_rule), prediction_strategy(prediction_strategy) {
}

void ForestModel::initCpp(uint mtry,
                     uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
                     std::string load_forest_filename, uint min_node_size,
                     std::string split_select_weights_file, std::vector<std::string>& always_split_variable_names,
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

  // Call other init function
  init(mtry, num_trees, seed, num_threads,
       min_node_size, prediction_mode, sample_with_replacement,
       memory_saving_splitting, sample_fraction);

  // TODO: Read 2d weights for tree-wise split select weights
  // Load split select weights from file


  tree_model = new TreeModel(relabeling_strategy,
                             splitting_rule,
                             prediction_strategy,
                             dependent_varID,
                             mtry,
                             min_node_size,
                             &deterministic_varIDs,
                             &split_select_varIDs,
                             &no_split_variables);

}

void ForestModel::init(uint mtry,
                       uint num_trees, uint seed, uint num_threads,
                       uint min_node_size,
                       bool prediction_mode, bool sample_with_replacement,
                       bool memory_saving_splitting,
                       double sample_fraction) {
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

  // Init split select weights
  split_select_weights.push_back(std::vector<double>());
}


void ForestModel::writeOutput(Data* prediction_data, std::vector<std::vector<double>> predictions) {
  *verbose_out << std::endl;
  *verbose_out << "Dependent variable name:           " << prediction_data->getVariableNames()[dependent_varID] << std::endl;
  *verbose_out << "Dependent variable ID:             " << dependent_varID << std::endl;
  *verbose_out << "Number of trees:                   " << num_trees << std::endl;
  *verbose_out << "Mtry:                              " << mtry << std::endl;
  *verbose_out << "Target node size:                  " << min_node_size << std::endl;
  *verbose_out << "Seed:                              " << seed << std::endl;
  *verbose_out << "Number of threads:                 " << num_threads << std::endl;
  *verbose_out << std::endl;

  if (prediction_mode) {
    writePredictionFile(prediction_data, predictions);
  } else {
    if (!split_select_weights.empty() & !split_select_weights[0].empty()) {
      *verbose_out
          << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
          << std::endl;
    }

    writeConfusionFile(prediction_data, predictions);
  }
}

Forest* ForestModel::train(Data* data) {
  size_t num_samples = data->getNumRows();
  size_t num_variables = data->getNumCols();
  size_t num_independent_variables = num_variables - no_split_variables.size();

  if (!split_select_weights_file.empty()) {
    std::vector<std::vector<double>> split_select_weights;
    split_select_weights.resize(1);
    loadDoubleVectorFromFile(split_select_weights[0], split_select_weights_file);
    if (split_select_weights[0].size() != num_variables - 1) {
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

  auto observations = new std::unordered_map<std::string, std::vector<double>>();
  for (auto it : observables) {
    std::string name = it.first;
    size_t index = it.second;

    for (int row = 0; row < data->getNumRows(); row++) {
      (*observations)[name].push_back(data->get(row, index));
    }
  }


  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);

  progress = 0;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::vector<std::future<std::vector<Tree*>>> futures;
  futures.reserve(num_threads);

  std::vector<Tree*> trees;
  trees.reserve(num_trees);

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::vector<Tree*>> promise;
    threads.push_back(std::thread(&ForestModel::growTreesInThread,
                                  this,
                                  i,
                                  data,
                                  observations,
                                  std::move(promise)));
    futures.push_back(promise.get_future());
  }
  showProgress("Growing trees..");
  for (auto &future : futures) {
    future.wait();
    std::vector<Tree*> thread_trees;
    trees.insert(trees.end(), thread_trees.begin(), thread_trees.end());
  }

  return new Forest(trees, data, observables);
}

std::vector<std::vector<double>> ForestModel::predict(Forest* forest, Data* prediction_data) {

  std::unordered_map<Tree*, std::vector<size_t>> terminal_node_IDs_by_tree;

  progress = 0;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::vector<std::future<
      std::unordered_map<Tree*, std::vector<size_t>>>> futures;

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::unordered_map<Tree *, std::vector<size_t>>> promise;
    threads.push_back(std::thread(&ForestModel::predictTreesInThread,
                                  this,
                                  i,
                                  forest,
                                  prediction_data,
                                  false,
                                  std::move(promise)));
    futures.push_back(promise.get_future());
  }

  showProgress("Predicting..");
  for (auto &future : futures) {
    future.wait();
    std::unordered_map<Tree*, std::vector<size_t>> terminal_nodeIDs = future.get();
    terminal_node_IDs_by_tree.insert(terminal_nodeIDs.begin(), terminal_nodeIDs.end());
  }

  return predictInternal(forest, prediction_data, terminal_node_IDs_by_tree);
}

void ForestModel::computePredictionError(Forest* forest,
                                         Data* prediction_data) {
// Predict trees in multiple threads
  std::unordered_map<Tree*, std::vector<size_t>> terminal_node_IDs_by_tree;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::vector<std::future<
      std::unordered_map<Tree*, std::vector<size_t>>>> futures;

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::unordered_map<Tree *, std::vector<size_t>>> promise;
    threads.push_back(std::thread(&ForestModel::predictTreesInThread,
                                  this,
                                  i,
                                  forest,
                                  prediction_data,
                                  false,
                                  std::move(promise)));
    futures.push_back(promise.get_future());
  }
  showProgress("Computing prediction error..");
  for (auto &future : futures) {
    future.wait();
    std::unordered_map<Tree*, std::vector<size_t>> terminal_nodeIDs = future.get();
    terminal_node_IDs_by_tree.insert(terminal_nodeIDs.begin(), terminal_nodeIDs.end());
  }
  computePredictionErrorInternal(forest, prediction_data, terminal_node_IDs_by_tree);
}

void ForestModel::growTreesInThread(uint thread_idx,
                                    Data* data,
                                    std::unordered_map<std::string, std::vector<double>>* observations,
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

      // Get split select weights for tree
      std::vector<double>* tree_split_select_weights;
      if (split_select_weights.size() > 1) {
        tree_split_select_weights = &split_select_weights[i];
      } else {
        tree_split_select_weights = &split_select_weights[0];
      }

      BootstrapSampler* bootstrap_sampler = new BootstrapSampler(data->getNumRows(),
                                                                 tree_seed,
                                                                 sample_with_replacement,
                                                                 sample_fraction,
                                                                 keep_inbag,
                                                                 &case_weights);
      trees.push_back(tree_model->train(data,
                                        observations,
                                        bootstrap_sampler,
                                        tree_split_select_weights));

      // Increase progress by 1 tree
      ++progress;
      condition_variable.notify_one();
    }
  }
  promise.set_value(trees);
}

void ForestModel::predictTreesInThread(uint thread_idx,
                                       Forest *forest,
                                       const Data *prediction_data,
                                       bool oob_prediction,
                                       std::promise<std::unordered_map<Tree*, std::vector<size_t>>> promise) {
  std::unordered_map<Tree *, std::vector<size_t>> terminal_nodeIDs_by_tree;

  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      Tree *tree = forest->get_trees()[i];
      std::vector<size_t> terminal_nodeIDs = tree->predict(prediction_data, oob_prediction);
      terminal_nodeIDs_by_tree[tree] = terminal_nodeIDs;

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
  promise.set_value(terminal_nodeIDs_by_tree);
}

void ForestModel::setSplitWeightVector(std::vector<std::vector<double>> &split_select_weights,
                                       size_t num_independent_variables) {

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

void ForestModel::setAlwaysSplitVariables(Data* data,
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

void ForestModel::showProgress(std::string operation) {
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

std::vector<std::vector<double>> ForestModel::predictInternal(Forest* forest,
                                                              Data* prediction_data,
                                                              std::unordered_map<Tree*, std::vector<size_t>> terminal_node_IDs_by_tree) {
  std::vector<std::vector<double>> predictions;
  predictions.reserve(prediction_data->getNumRows());

  for (size_t sampleID = 0; sampleID < prediction_data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < forest->get_trees().size(); ++tree_idx) {
      Tree* tree = forest->get_trees()[tree_idx];
      std::vector<size_t> terminal_node_IDs = terminal_node_IDs_by_tree[tree];
      addSampleWeights(sampleID, tree, weights_by_sampleID, terminal_node_IDs);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> prediction = prediction_strategy->predict(weights_by_sampleID, forest->get_original_observations());
    predictions.push_back(prediction);
  }
  return predictions;
}

void ForestModel::addSampleWeights(size_t test_sample_idx,
                                   Tree *tree,
                                   std::unordered_map<size_t, double> &weights_by_sampleID,
                                   std::vector<size_t> terminal_node_IDs) {
  size_t nodeID = terminal_node_IDs[test_sample_idx];
  std::vector<size_t> sampleIDs = tree->get_sampleIDs()[nodeID];
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestModel::normalizeSampleWeights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}

void ForestModel::computePredictionErrorInternal(Forest* forest,
                                                 Data* prediction_data,
                                                 std::unordered_map<Tree*, std::vector<size_t>> terminal_node_IDs_by_tree) {
  std::vector<std::vector<double>> predictions;
  predictions.reserve(prediction_data->getNumRows());

  std::unordered_multimap<size_t, size_t> trees_by_oob_samples;
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sampleID : forest->get_trees()[tree_idx]->getOobSampleIDs()) {
      trees_by_oob_samples.insert(std::pair<size_t, size_t>(sampleID, tree_idx));
    }
  }

  for (size_t sampleID = 0; sampleID < prediction_data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;

    auto tree_range = trees_by_oob_samples.equal_range(sampleID);
    if (tree_range.first == tree_range.second) {
      std::vector<double> temp{0.0};
      predictions.push_back(temp);
      continue;
    }

    // Calculate the weights of neighboring samples.
    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      Tree *tree = forest->get_trees()[it->second];

      // hackhackhack
      std::vector<size_t> oob_sampleIDs = tree->getOobSampleIDs();
      size_t sample_idx = (size_t) (std::find(oob_sampleIDs.begin(), oob_sampleIDs.end(), sampleID) -
                                    oob_sampleIDs.begin());

      addSampleWeights(sample_idx, tree, weights_by_sampleID, terminal_node_IDs_by_tree[tree]);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> prediction = prediction_strategy->predict(weights_by_sampleID, forest->get_original_observations());
    predictions.push_back(prediction);
  }
}

void ForestModel::writeConfusionFile(Data* prediction_data, std::vector<std::vector<double>> predictions) {

// Open confusion file for writing
  std::string filename = "ranger.confusion";
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

void ForestModel::writePredictionFile(Data* prediction_data, std::vector<std::vector<double>> predictions) {
  std::string filename = "ranger.prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Prediction X1" << std::endl;
  for (size_t i = 0; i < predictions.size(); ++i) {
    for (size_t j = 0; j < predictions[i].size(); ++j) {
      outfile << predictions[i][j] << " " << prediction_data->get(i, 0) << " ";
    }
    outfile << std::endl;
  }

  *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}
