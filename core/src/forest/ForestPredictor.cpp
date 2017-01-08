#include "ForestPredictor.h"
#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <thread>
#include <future>
#include "utility.h"
#include "ForestPredictor.h"

ForestPredictor::ForestPredictor(uint num_threads,
                                 std::shared_ptr<PredictionStrategy> prediction_strategy) :
    prediction_strategy(prediction_strategy) {
  if (num_threads == DEFAULT_NUM_THREADS) {
    this->num_threads = std::thread::hardware_concurrency();
  } else {
    this->num_threads = num_threads;
  }
}

std::vector<std::vector<double>> ForestPredictor::predict(const Forest& forest, Data* prediction_data) {
  std::unordered_map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree;

  size_t num_trees = forest.get_trees().size();
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::vector<std::future<
  std::unordered_map<size_t, std::vector<size_t>>>> futures;

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::unordered_map<size_t, std::vector<size_t>>> promise;
    std::future<std::unordered_map<size_t, std::vector<size_t>>> future = promise.get_future();
    threads.push_back(std::thread(&ForestPredictor::predictTreesInThread,
                                  this,
                                  i,
                                  forest,
                                  prediction_data,
                                  false,
                                  std::move(promise)));
    futures.push_back(std::move(future));
  }

  for (auto &future : futures) {
    future.wait();
    std::unordered_map<size_t, std::vector<size_t>> terminal_nodeIDs = future.get();
    terminal_node_IDs_by_tree.insert(terminal_nodeIDs.begin(), terminal_nodeIDs.end());
  }

  for (auto& thread : threads) {
    thread.join();
  }

  return predictInternal(forest, prediction_data, terminal_node_IDs_by_tree);
}

void ForestPredictor::computePredictionError(const Forest& forest,
                                             Data* prediction_data) {
// Predict trees in multiple threads
  std::unordered_map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::vector<std::future<
  std::unordered_map<size_t, std::vector<size_t>>>> futures;

  for (uint i = 0; i < num_threads; ++i) {
    std::promise<std::unordered_map<size_t, std::vector<size_t>>> promise;
    std::future<std::unordered_map<size_t, std::vector<size_t>>> future = promise.get_future();
    threads.push_back(std::thread(&ForestPredictor::predictTreesInThread,
                                  this,
                                  i,
                                  forest,
                                  prediction_data,
                                  false,
                                  std::move(promise)));
    futures.push_back(std::move(future));
  }

  for (auto &future : futures) {
    future.wait();
    std::unordered_map<size_t, std::vector<size_t>> terminal_nodeIDs = future.get();
    terminal_node_IDs_by_tree.insert(terminal_nodeIDs.begin(), terminal_nodeIDs.end());
  }

  for (auto& thread: threads) {
    thread.join();
  }

  computePredictionErrorInternal(forest, prediction_data, terminal_node_IDs_by_tree);
}

void ForestPredictor::predictTreesInThread(uint thread_idx,
                                           const Forest& forest,
                                           const Data *prediction_data,
                                           bool oob_prediction,
                                           std::promise<std::unordered_map<size_t, std::vector<size_t>>> promise) {
  std::unordered_map<size_t, std::vector<size_t>> terminal_nodeIDs_by_tree;

  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      std::shared_ptr<Tree> tree = forest.get_trees()[i];
      std::vector<size_t> terminal_nodeIDs = tree->predict(prediction_data, oob_prediction);
      terminal_nodeIDs_by_tree[i] = terminal_nodeIDs;
    }
  }
  promise.set_value(terminal_nodeIDs_by_tree);
}

std::vector<std::vector<double>> ForestPredictor::predictInternal(const Forest& forest,
      Data* prediction_data,
      const std::unordered_map<size_t, std::vector<size_t>>& terminal_node_IDs_by_tree) {
  std::vector<std::vector<double>> predictions;
  predictions.reserve(prediction_data->getNumRows());

  for (size_t sampleID = 0; sampleID < prediction_data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < forest.get_trees().size(); ++tree_idx) {
      std::shared_ptr<Tree> tree = forest.get_trees()[tree_idx];
      std::vector<size_t> terminal_node_IDs = terminal_node_IDs_by_tree.at(tree_idx);
      addSampleWeights(sampleID, tree, weights_by_sampleID, terminal_node_IDs);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> prediction = prediction_strategy->predict(weights_by_sampleID,
                                                                  forest.get_observations());
    predictions.push_back(prediction);
  }
  return predictions;
}

void ForestPredictor::addSampleWeights(size_t test_sample_idx,
                                       std::shared_ptr<Tree> tree,
                                       std::unordered_map<size_t, double> &weights_by_sampleID,
                                       std::vector<size_t> terminal_node_IDs) {
  size_t nodeID = terminal_node_IDs[test_sample_idx];
  std::vector<size_t> sampleIDs = tree->get_sampleIDs()[nodeID];
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestPredictor::normalizeSampleWeights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}

void ForestPredictor::computePredictionErrorInternal(const Forest& forest,
    Data* prediction_data,
    const std::unordered_map<size_t, std::vector<size_t>>& terminal_node_IDs_by_tree) {
  std::vector<std::vector<double>> predictions;
  predictions.reserve(prediction_data->getNumRows());

  size_t num_trees = forest.get_trees().size();
  std::unordered_multimap<size_t, size_t> trees_by_oob_samples;
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sampleID : forest.get_trees()[tree_idx]->getOobSampleIDs()) {
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
      std::shared_ptr<Tree> tree = forest.get_trees()[it->second];

      // hackhackhack
      std::vector<size_t> oob_sampleIDs = tree->getOobSampleIDs();
      size_t sample_idx = (size_t) (std::find(oob_sampleIDs.begin(), oob_sampleIDs.end(), sampleID) -
                                    oob_sampleIDs.begin());

      addSampleWeights(sample_idx, tree, weights_by_sampleID, terminal_node_IDs_by_tree.at(it->second));
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> prediction = prediction_strategy->predict(weights_by_sampleID,
                                                                  forest.get_observations());
    predictions.push_back(prediction);
  }
}

void ForestPredictor::writeConfusionFile(Data* prediction_data, std::vector<std::vector<double>> predictions) {

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
}

void ForestPredictor::writePredictionFile(Data* prediction_data, std::vector<std::vector<double>> predictions) {
  std::string filename = "gradientforest.prediction";
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
}
