#include <algorithm>
#include <stdexcept>
#include <string>
#include <Tree/QuantileRelabelingStrategy.h>
#include <set>
#include <Tree/ProbabilitySplittingRule.h>

#include "utility.h"
#include "ForestQuantile.h"
#include "PredictionStrategy.h"
#include "QuantilePredictionStrategy.h"

ForestQuantile::ForestQuantile(std::unordered_map<std::string, size_t> observables,
                               RelabelingStrategy *relabeling_strategy,
                               SplittingRule *splitting_rule,
                               PredictionStrategy* prediction_strategy):
    Forest(observables, relabeling_strategy, splitting_rule, prediction_strategy) {}

ForestQuantile::~ForestQuantile() {}

void ForestQuantile::initInternal(std::string status_variable_name) {

  // If mtry not set, use number of independent variables / 3.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }
}

void ForestQuantile::predictInternal() {
  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
      Tree* tree = trees[tree_idx];
      addSampleWeights(sampleID, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> quantile_cutoffs = prediction_strategy->predict(weights_by_sampleID, original_observations);
    predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::addSampleWeights(size_t test_sample_idx,
                                      Tree* tree,
                                      std::unordered_map<size_t, double> &weights_by_sampleID) {
  std::vector<size_t> sampleIDs = tree->get_neighboring_samples(test_sample_idx);
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestQuantile::normalizeSampleWeights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}

void ForestQuantile::computePredictionErrorInternal() {
  std::unordered_multimap<size_t, size_t> trees_by_oob_samples;
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sampleID : trees[tree_idx]->getOobSampleIDs()) {
      trees_by_oob_samples.insert(std::pair<size_t, size_t>(sampleID, tree_idx));
    }
  }

  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;

    auto tree_range = trees_by_oob_samples.equal_range(sampleID);
    if (tree_range.first == tree_range.second) {
      predictions.push_back(std::vector<double>(1, -1.0));
      continue;
    }

    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      Tree* tree = trees[it->second];

      // hackhackhack
      std::vector<size_t> oob_sampleIDs = tree->getOobSampleIDs();
      size_t sample_idx = (size_t) (std::find(oob_sampleIDs.begin(), oob_sampleIDs.end(), sampleID) - oob_sampleIDs.begin());

      addSampleWeights(sample_idx, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<std::pair<size_t, double>> sampleIDs_and_values;
    for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
      size_t sampleID = it->first;
      sampleIDs_and_values.push_back(std::pair<size_t, double>(
          sampleID, data->get(sampleID, dependent_varID)));
    }

    std::vector<double> quantile_cutoffs = prediction_strategy->predict(weights_by_sampleID, original_observations);
    predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Quantile" << std::endl;
}

void ForestQuantile::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = "ranger.confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  for (auto &quantile_prediction : predictions) {
    for (auto it = quantile_prediction.begin(); it != quantile_prediction.end(); ++it) {
      outfile << *it << " ";
    }
    outfile << std::endl;
  }

  outfile.close();
  *verbose_out << "Saved prediction error to file " << filename << "." << std::endl;
}

void ForestQuantile::writePredictionFile() {

// Open prediction file for writing
  std::string filename = "ranger.prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Predictions: " << std::endl;
  for (size_t i = 0; i < predictions.size(); ++i) {
    for (size_t j = 0; j < predictions[i].size(); ++j) {
      outfile << predictions[i][j] << " ";
    }
    outfile << std::endl;
  }

  *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}
