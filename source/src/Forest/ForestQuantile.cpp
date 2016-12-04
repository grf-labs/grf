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

ForestQuantile::ForestQuantile(std::vector<double>* quantiles,
                               RelabelingStrategy *relabeling_strategy,
                               SplittingRule *splitting_rule):
    Forest(relabeling_strategy, splitting_rule),
    quantiles(quantiles),
    original_responses(new std::vector<double>()) {}

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

  if (!prediction_mode) {
    original_responses = new std::vector<double>();
    for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
      original_responses->push_back(data->get(sampleID, dependent_varID));
    }
  }
}

void ForestQuantile::predictInternal() {
  PredictionStrategy* predictionStrategy = new QuantilePredictionStrategy(quantiles, original_responses);
  
  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
      TreeFactory* tree = trees[tree_idx];
      addSampleWeights(sampleID, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> quantile_cutoffs = predictionStrategy->predict(weights_by_sampleID);
    predictions.push_back(quantile_cutoffs);
  }
}

void ForestQuantile::addSampleWeights(size_t test_sample_idx,
                                      TreeFactory* tree,
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
  PredictionStrategy* predictionStrategy = new QuantilePredictionStrategy(quantiles, original_responses);

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
      predictions.push_back(std::vector<double>(quantiles->size(), 0.0));
      continue;
    }

    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      TreeFactory* tree = trees[it->second];

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

    std::vector<double> quantile_cutoffs = predictionStrategy->predict(weights_by_sampleID);
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

void ForestQuantile::saveToFileInternal(std::ofstream& outfile) {

// Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

// Write treetype
  TreeType treetype = TREE_QUANTILE;
  outfile.write((char*) &treetype, sizeof(treetype));
  saveVector1D(*quantiles, outfile);
  saveVector1D(*original_responses, outfile);
}

void ForestQuantile::loadFromFileInternal(std::ifstream& infile) {

// Read number of variables
  size_t num_variables_saved;
  infile.read((char *) &num_variables_saved, sizeof(num_variables_saved));

// Read treetype
  TreeType treetype;
  infile.read((char *) &treetype, sizeof(treetype));
  if (treetype != TREE_QUANTILE) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a regression forest.");
  }

  readVector1D(*quantiles, infile);
  readVector1D(*original_responses, infile);

  for (size_t i = 0; i < num_trees; ++i) {

    // Read data
    std::vector<std::vector<size_t>> child_nodeIDs;
    readVector2D(child_nodeIDs, infile);
    std::vector<size_t> split_varIDs;
    readVector1D(split_varIDs, infile);
    std::vector<double> split_values;
    readVector1D(split_values, infile);

    // If dependent variable not in test data, change variable IDs accordingly
    if (num_variables_saved > num_variables) {
      for (auto &varID : split_varIDs) {
        if (varID >= dependent_varID) {
          --varID;
        }
      }
    }

    std::vector<std::vector<size_t>> sampleIDs;
    readVector2D(sampleIDs, infile);

    // Create tree
    RelabelingStrategy* relabeling_strategy = new QuantileRelabelingStrategy(
        quantiles,
        dependent_varID);
    SplittingRule *splitting_rule = new ProbabilitySplittingRule(data, quantiles->size());

    TreeFactory *tree = new TreeFactory(child_nodeIDs, split_varIDs, split_values,
                                        sampleIDs, relabeling_strategy, splitting_rule);
    trees.push_back(tree);
  }
}