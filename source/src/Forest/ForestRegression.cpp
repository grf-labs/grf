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

#include <algorithm>
#include <stdexcept>
#include <string>

#include "utility.h"
#include "ForestRegression.h"
#include "TreeRegression.h"
#include "Data.h"

ForestRegression::ForestRegression() {
}

ForestRegression::~ForestRegression() {
}

void ForestRegression::loadForest(size_t dependent_varID, size_t num_trees,
    std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
    std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
    std::vector<bool>& is_ordered_variable) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->is_ordered_variable = is_ordered_variable;

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    Tree* tree = new TreeRegression(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
        &this->is_ordered_variable);
    trees.push_back(tree);
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestRegression::initInternal(std::string status_variable_name) {

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

void ForestRegression::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeRegression());
  }
}

void ForestRegression::predictInternal() {

  size_t num_prediction_samples = data->getNumRows();
  predictions.reserve(num_prediction_samples);

  // For all samples get tree predictions
  for (size_t sample_idx = 0; sample_idx < num_prediction_samples; ++sample_idx) {

    if (predict_all) {
      // Get all tree predictions
      std::vector<double> sample_predictions;
      sample_predictions.reserve(num_trees);
      for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
        double value = ((TreeRegression*) trees[tree_idx])->getPrediction(sample_idx);
        sample_predictions.push_back(value);
      }
      predictions.push_back(sample_predictions);
    } else {
      // Mean over trees
      double prediction_sum = 0;
      for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
        prediction_sum += ((TreeRegression*) trees[tree_idx])->getPrediction(sample_idx);
      }
      std::vector<double> temp;
      temp.push_back(prediction_sum / num_trees);
      predictions.push_back(temp);
    }
  }
}

void ForestRegression::computePredictionErrorInternal() {

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
      double value = ((TreeRegression*) trees[tree_idx])->getPrediction(sample_idx);

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

void ForestRegression::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Regression" << std::endl;
}

void ForestRegression::writeConfusionFile() {

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

void ForestRegression::writePredictionFile() {

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

void ForestRegression::saveToFileInternal(std::ofstream& outfile) {

// Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

// Write treetype
  TreeType treetype = TREE_REGRESSION;
  outfile.write((char*) &treetype, sizeof(treetype));
}

void ForestRegression::loadFromFileInternal(std::ifstream& infile) {

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
    Tree* tree = new TreeRegression(child_nodeIDs, split_varIDs, split_values, &is_ordered_variable);
    trees.push_back(tree);
  }
}

