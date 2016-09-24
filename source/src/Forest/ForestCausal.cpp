#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestCausal.h"
#include "TreeCausal.h"

ForestCausal::ForestCausal(): treatment_varID(0) {}

ForestCausal::~ForestCausal() {}

void ForestCausal::loadForest(size_t dependent_varID, size_t num_trees,
                              std::vector<std::vector<std::vector<size_t>>> &forest_child_nodeIDs,
                              std::vector<std::vector<size_t>> &forest_split_varIDs,
                              std::vector<std::vector<double>> &forest_split_values) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    Tree* tree = new TreeCausal(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
                                std::vector<std::vector<size_t>>(), treatment_varID);
    trees.push_back(tree);
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestCausal::initInternal(std::string status_variable_name) {
  if (!prediction_mode) {
    treatment_varID = data->getVariableID(status_variable_name);
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

void ForestCausal::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeCausal(treatment_varID));
  }
}

void ForestCausal::predictInternal() {
  size_t num_prediction_samples = data->getNumRows();
  predictions.reserve(num_prediction_samples);

  // For all samples get tree predictions
  for (size_t sample_idx = 0; sample_idx < num_prediction_samples; ++sample_idx) {

    if (predict_all) {
      throw std::runtime_error("Causal forests do not support predict_all.");
    } else {
      // Mean over trees
      double prediction_sum = 0;
      for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
        prediction_sum += ((TreeCausal*) trees[tree_idx])->getPrediction(sample_idx);
      }
      std::vector<double> temp;
      temp.push_back(prediction_sum / num_trees);
      predictions.push_back(temp);
    }
  }
}

void ForestCausal::computePredictionErrorInternal() {

// For each sample sum over trees where sample is OOB
  std::vector<size_t> samples_oob_count;
  predictions.reserve(num_samples);
  samples_oob_count.resize(num_samples, 0);
  for (size_t i = 0; i < num_samples; ++i) {
    std::vector<double> temp { 0 };
    predictions.push_back(temp);
  }
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sample_idx = 0; sample_idx < trees[tree_idx]->getNumSamplesOob(); ++sample_idx) {
      size_t sampleID = trees[tree_idx]->getOobSampleIDs()[sample_idx];
      double value = ((TreeCausal*) trees[tree_idx])->getPrediction(sample_idx);

      predictions[sampleID][0] += value;
      ++samples_oob_count[sampleID];
    }
  }

  for (size_t i = 0; i < predictions.size(); ++i) {
    if (samples_oob_count[i] > 0) {
      predictions[i][0] /= (double) samples_oob_count[i];
    }
  }
}

void ForestCausal::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Regression" << std::endl;
}

void ForestCausal::writeConfusionFile() {

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

void ForestCausal::writePredictionFile() {

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

void ForestCausal::saveToFileInternal(std::ofstream& outfile) {

// Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

// Write treetype
  TreeType treetype = TREE_CAUSAL;
  outfile.write((char*) &treetype, sizeof(treetype));

  outfile.write((char*) &treatment_varID, sizeof(treatment_varID));
}

void ForestCausal::loadFromFileInternal(std::ifstream& infile) {

// Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

// Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_CAUSAL) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a causal forest.");
  }

  infile.read((char*) &treatment_varID, sizeof(treatment_varID));

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
    Tree* tree = new TreeCausal(child_nodeIDs, split_varIDs, split_values,
    std::vector<std::vector<size_t>>(), treatment_varID);
    trees.push_back(tree);
  }
}