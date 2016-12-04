#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestInstrumental.h"
#include "InstrumentalRelabelingStrategy.h"
#include "RegressionSplittingRule.h"
#include "PredictionStrategy.h"
#include "InstrumentalPredictionStrategy.h"

ForestInstrumental::ForestInstrumental(RelabelingStrategy* relabeling_strategy,
                                       SplittingRule* splitting_rule,
                                       std::string instrument_variable_name) :
    Forest(relabeling_strategy, splitting_rule),
    treatment_varID(0),
    instrument_varID(0),
    instrument_variable_name(instrument_variable_name),
    original_responses(0) {}

ForestInstrumental::~ForestInstrumental() {}

void ForestInstrumental::initInternal(std::string status_variable_name) {
  if (!prediction_mode) {
    treatment_varID = data->getVariableID(status_variable_name);
    instrument_varID = data->getVariableID(instrument_variable_name);

    no_split_variables.push_back(treatment_varID);
    no_split_variables.push_back(instrument_varID);
  }

  // If mtry not set, use the square root of the number of possible split variables.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - no_split_variables.size()));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }

  if (!prediction_mode) {
    original_responses = new std::unordered_map<size_t, std::vector<double>>;

    for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
      (*original_responses)[dependent_varID].push_back(data->get(sampleID, dependent_varID));
    }

    for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
      (*original_responses)[dependent_varID + 1].push_back(data->get(sampleID, treatment_varID));
    }

    for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
      (*original_responses)[dependent_varID + 2].push_back(data->get(sampleID, instrument_varID));
    }
  }
}

void ForestInstrumental::predictInternal() {
  PredictionStrategy* predictionStrategy = new InstrumentalPredictionStrategy(instrument_varID,
                                                                              treatment_varID,
                                                                              dependent_varID,
                                                                              original_responses);

  predictions.reserve(num_samples);

  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
    std::unordered_map<size_t, double> weights_by_sampleID;
    for (size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
      TreeFactory* tree = trees[tree_idx];
      addSampleWeights(sampleID, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> prediction = predictionStrategy->predict(weights_by_sampleID);
    predictions.push_back(prediction);
  }
}

void ForestInstrumental::addSampleWeights(size_t test_sample_idx,
                                          TreeFactory *tree,
                                          std::unordered_map<size_t, double> &weights_by_sampleID) {
  std::vector<size_t> sampleIDs = tree->get_neighboring_samples(test_sample_idx);
  double sample_weight = 1.0 / sampleIDs.size();

  for (auto &sampleID : sampleIDs) {
    weights_by_sampleID[sampleID] += sample_weight;
  }
}

void ForestInstrumental::normalizeSampleWeights(std::unordered_map<size_t, double>& weights_by_sampleID) {
  double total_weight = 0.0;
  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    total_weight += it->second;
  }

  for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
    it->second /= total_weight;
  }
}

void ForestInstrumental::computePredictionErrorInternal() {
  PredictionStrategy* predictionStrategy = new InstrumentalPredictionStrategy(instrument_varID,
                                                                              treatment_varID,
                                                                              dependent_varID,
                                                                              original_responses);
  predictions.reserve(num_samples);

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
      std::vector<double> temp { 0.0 };
      predictions.push_back(temp);
      continue;
    }

    // Calculate the weights of neighboring samples.
    for (auto it = tree_range.first; it != tree_range.second; ++it) {
      TreeFactory* tree = trees[it->second];

      // hackhackhack
      std::vector<size_t> oob_sampleIDs = tree->getOobSampleIDs();
      size_t sample_idx = (size_t) (std::find(oob_sampleIDs.begin(), oob_sampleIDs.end(), sampleID) - oob_sampleIDs.begin());

      addSampleWeights(sample_idx, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    std::vector<double> prediction = predictionStrategy->predict(weights_by_sampleID);
    predictions.push_back(prediction);
  }
}

void ForestInstrumental::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Regression" << std::endl;
}

void ForestInstrumental::writeConfusionFile() {

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
      outfile << predictions[i][j] << " " << data->get(i, 0) << " ";
    }
    outfile << std::endl;
  }

  outfile.close();
  *verbose_out << "Saved prediction error to file " << filename << "." << std::endl;
}

void ForestInstrumental::writePredictionFile() {

// Open prediction file for writing
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
      outfile << predictions[i][j] << " " << data->get(i, 0) << " ";
    }
    outfile << std::endl;
  }

  *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}

void ForestInstrumental::saveToFileInternal(std::ofstream& outfile) {

// Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

// Write treetype
  TreeType treetype = TREE_INSTRUMENTAL;
  outfile.write((char*) &treetype, sizeof(treetype));

  outfile.write((char*) &treatment_varID, sizeof(treatment_varID));
  outfile.write((char*) &instrument_varID, sizeof(instrument_varID));

  saveVector1D((*original_responses)[dependent_varID], outfile);
  saveVector1D((*original_responses)[treatment_varID], outfile);
  saveVector1D((*original_responses)[instrument_varID], outfile);
}

void ForestInstrumental::loadFromFileInternal(std::ifstream& infile) {

// Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

// Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_INSTRUMENTAL) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a instrumental forest.");
  }

  infile.read((char*) &treatment_varID, sizeof(treatment_varID));
  infile.read((char*) &instrument_varID, sizeof(instrument_varID));

  original_responses = new std::unordered_map<size_t, std::vector<double>>();
  readVector1D((*original_responses)[dependent_varID], infile);
  readVector1D((*original_responses)[treatment_varID], infile);
  readVector1D((*original_responses)[instrument_varID], infile);

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
      for (auto& varID : split_varIDs) {
        if (varID >= dependent_varID) {
          --varID;
        }
      }
    }

    std::vector<std::vector<size_t>> sampleIDs;
    readVector2D(sampleIDs, infile);

    RelabelingStrategy* relabeling_strategy = new InstrumentalRelabelingStrategy(
        dependent_varID, treatment_varID, instrument_varID);
    SplittingRule* splitting_rule = new RegressionSplittingRule(data);
    TreeFactory *tree = new TreeFactory(child_nodeIDs, split_varIDs, split_values,
                                      sampleIDs, relabeling_strategy, splitting_rule);
    trees.push_back(tree);
  }
}