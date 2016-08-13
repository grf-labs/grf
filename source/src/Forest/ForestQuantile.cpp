#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestQuantile.h"
#include "TreeQuantile.h"
#include "Data.h"

ForestQuantile::ForestQuantile(std::vector<double>* quantiles): quantiles(quantiles) {}

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

  // Sort data if memory saving mode
  if (!memory_saving_splitting) {
    data->sort();
  }
}

void ForestQuantile::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeQuantile(quantiles));
  }
}

void ForestQuantile::predictInternal() {

  size_t num_prediction_samples = data->getNumRows();
  predictions.reserve(num_prediction_samples);

  // For all samples get tree predictions
  for (size_t sample_idx = 0; sample_idx < num_prediction_samples; ++sample_idx) {
    if (predict_all) {
      throw std::runtime_error("Quantile forests do not support 'predict_all'.");
    }

    double prediction_sum = 0;
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      prediction_sum += ((TreeQuantile*) trees[tree_idx])->getPrediction(sample_idx);
    }
    std::vector<double> temp;
    temp.push_back(prediction_sum / num_trees);
    predictions.push_back(temp);
  }
}

void ForestQuantile::computePredictionErrorInternal() {

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
      double value = ((TreeQuantile*) trees[tree_idx])->getPrediction(sample_idx);

      predictions[sampleID][0] += value;
      ++samples_oob_count[sampleID];
    }
  }

// MSE with predictions and true data
//oob_anytree_sampleIDs.reserve(predictions.size());
  for (size_t i = 0; i < predictions.size(); ++i) {
    if (samples_oob_count[i] > 0) {
      //oob_anytree_sampleIDs.push_back(i);
      predictions[i][0] /= (double) samples_oob_count[i];
      double predicted_value = predictions[i][0];
      double real_value = data->get(i, dependent_varID);
      overall_prediction_error += (predicted_value - real_value) * (predicted_value - real_value);
    }
  }

  overall_prediction_error /= (double) predictions.size();
}

void ForestQuantile::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Regression" << std::endl;
}

void ForestQuantile::writeConfusionFile() {

// Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

// Write confusion to file
  outfile << "Overall OOB prediction error (MSE): " << overall_prediction_error << std::endl;

  outfile.close();
  *verbose_out << "Saved prediction error to file " << filename << "." << std::endl;
}

void ForestQuantile::writePredictionFile() {

// Open prediction file for writing
  std::string filename = output_prefix + ".prediction";
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
  TreeType treetype = TREE_REGRESSION;
  outfile.write((char*) &treetype, sizeof(treetype));
}

void ForestQuantile::loadFromFileInternal(std::ifstream& infile) {

// Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

// Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_REGRESSION) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a regression forest.");
  }

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
    Tree* tree = new TreeQuantile(child_nodeIDs, split_varIDs, split_values, &is_ordered_variable);
    trees.push_back(tree);
  }
}