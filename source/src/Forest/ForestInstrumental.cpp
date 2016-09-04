#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestInstrumental.h"

ForestInstrumental::ForestInstrumental(std::string instrument_variable_name) :
    treatment_varID(0),
    instrument_varID(0),
    instrument_variable_name(instrument_variable_name),
    original_responses(0) {}

ForestInstrumental::~ForestInstrumental() {}

void ForestInstrumental::loadForest(size_t dependent_varID, size_t num_trees,
                                    std::vector<std::vector<std::vector<size_t>>> &forest_child_nodeIDs,
                                    std::vector<std::vector<size_t>> &forest_split_varIDs,
                                    std::vector<std::vector<double>> &forest_split_values,
                                    std::vector<bool> &is_ordered_variable,
                                    std::vector<double> *original_responses) {
  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->is_ordered_variable = is_ordered_variable;
  this->original_responses = new std::vector<double>(*original_responses);

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    Tree *tree = new TreeInstrumental(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
                                      &this->is_ordered_variable, treatment_varID, instrument_varID);
    trees.push_back(tree);
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestInstrumental::initInternal(std::string status_variable_name) {
  if (!prediction_mode) {
    treatment_varID = data->getVariableID(status_variable_name);
    instrument_varID = data->getVariableID(instrument_variable_name);
  }

  // If mtry not set, use number of independent variables / 3.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }

  // Sort data if memory saving mode
  if (!memory_saving_splitting) {
    data->sort();
  }
}

void ForestInstrumental::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeInstrumental(treatment_varID, instrument_varID));
  }
}

void ForestInstrumental::predictInternal() {
  if (predict_all) {
    throw std::runtime_error("Instrumental forests do not support 'predict_all'.");
  }

//  for (size_t sampleID = 0; sampleID < data->getNumRows(); ++sampleID) {
//    std::unordered_map<size_t, double> weights_by_sampleID;
//    for (size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
//      TreeInstrumental* tree = (TreeInstrumental *) trees[tree_idx];
//      addSampleWeights(sampleID, tree, weights_by_sampleID);
//    }
//
//    normalizeSampleWeights(weights_by_sampleID);
//
//    double average_instrument = 0.0;
//    double average_treatment = 0.0;
//    double average_response = 0.0;
//
//    std::vector<std::pair<size_t, double>> sampleIDs_and_values;
//    for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
//      size_t sampleID = it->first;
//      sampleIDs_and_values.push_back(std::pair<size_t, double>(
//          sampleID, original_responses->at(sampleID)));
//    }
//  }
}

void ForestInstrumental::addSampleWeights(size_t test_sample_idx,
                                      TreeInstrumental* tree,
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
  if (predict_all) {
    throw std::runtime_error("Instrumental forests do not support 'predict_all'.");
  }

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
      TreeInstrumental* tree = (TreeInstrumental*) trees[it->second];

      // hackhackhack
      std::vector<size_t> oob_sampleIDs = tree->getOobSampleIDs();
      size_t sample_idx = (size_t) (std::find(oob_sampleIDs.begin(), oob_sampleIDs.end(), sampleID) - oob_sampleIDs.begin());

      addSampleWeights(sample_idx, tree, weights_by_sampleID);
    }

    normalizeSampleWeights(weights_by_sampleID);

    // Compute the relevant averages.
    double average_instrument = 0.0;
    double average_treatment = 0.0;
    double average_response = 0.0;

    for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
      size_t neighborID = it->first;
      double weight = it->second;

      average_instrument += weight * data->get(neighborID, instrument_varID);
      average_treatment += weight * data->get(neighborID, treatment_varID);
      average_response += weight * data->get(neighborID, dependent_varID);
    }

    // Finally, calculate the prediction.
    double numerator = 0.0;
    double denominator = 0.0;
    for (auto it = weights_by_sampleID.begin(); it != weights_by_sampleID.end(); ++it) {
      size_t neighborID = it->first;
      double weight = it->second;

      double instrument = data->get(neighborID, instrument_varID);
      double treatment = data->get(neighborID, treatment_varID);
      double response = weight * data->get(neighborID, dependent_varID);

      numerator += weight * (instrument - average_instrument) * (response - average_response);
      denominator += weight * (instrument - average_instrument) * (treatment - average_treatment);
    }

    std::vector<double> temp { numerator / denominator };
    predictions.push_back(temp);
  }
}

void ForestInstrumental::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Regression" << std::endl;
}

void ForestInstrumental::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
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
  std::string filename = output_prefix + ".prediction";
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
  TreeType treetype = TREE_CAUSAL;
  outfile.write((char*) &treetype, sizeof(treetype));

  outfile.write((char*) &treatment_varID, sizeof(treatment_varID));
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

    // Create tree
    Tree *tree = new TreeInstrumental(child_nodeIDs, split_varIDs, split_values,
                                      &is_ordered_variable, treatment_varID, instrument_varID);
    trees.push_back(tree);
  }
}