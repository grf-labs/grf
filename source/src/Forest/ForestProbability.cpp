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

#include "utility.h"
#include "ForestProbability.h"
#include "TreeProbability.h"
#include "Data.h"

ForestProbability::ForestProbability() {
}

ForestProbability::~ForestProbability() {
}

void ForestProbability::loadForest(size_t dependent_varID, size_t num_trees,
    std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
    std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
    std::vector<double>& class_values, std::vector<std::vector<std::vector<double>>>& forest_terminal_class_counts) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->class_values = class_values;

  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    Tree* tree = new TreeProbability(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i],
    &class_values, &response_classIDs, forest_terminal_class_counts[i]);
    trees.push_back(tree);
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestProbability::initInternal(std::string status_variable_name) {

  // If mtry not set, use floored square root of number of independent variables.
  if (mtry == 0) {
    int temp = sqrt(num_variables - 1);
    mtry = std::max(1, temp);
  }

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_PROBABILITY;
  }

  // Create class_values and response_classIDs
  if (!prediction_mode) {
    for (size_t i = 0; i < num_samples; ++i) {
      double value = data->get(i, dependent_varID);

      // If classID is already in class_values, use ID. Else create a new one.
      uint classID = find(class_values.begin(), class_values.end(), value) - class_values.begin();
      if (classID == class_values.size()) {
        class_values.push_back(value);
      }
      response_classIDs.push_back(classID);
    }
  }
}

void ForestProbability::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(new TreeProbability(&class_values, &response_classIDs));
  }
}

void ForestProbability::predictInternal() {

  // First dim samples, second dim classes
  size_t num_prediction_samples = trees[0]->getPredictions().size();
  predictions.resize(num_prediction_samples);
  for (size_t i = 0; i < num_prediction_samples; ++i) {
    predictions[i].resize(class_values.size(), 0);
  }

  // For all samples average proportions of trees
  for (size_t sample_idx = 0; sample_idx < num_prediction_samples; ++sample_idx) {

    // For each sample compute proportions in each tree and average over trees
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      std::vector<double> counts = trees[tree_idx]->getPredictions()[sample_idx];
      double sum = 0;
      for (size_t class_idx = 0; class_idx < counts.size(); ++class_idx) {
        sum += counts[class_idx];
      }
      for (size_t class_idx = 0; class_idx < counts.size(); ++class_idx) {
        predictions[sample_idx][class_idx] += counts[class_idx] / sum / num_trees;
      }
    }
  }

}

void ForestProbability::computePredictionErrorInternal() {

  // Class counts for samples
  std::vector<std::vector<size_t>> class_counts;
  class_counts.resize(num_samples, std::vector<size_t>());
  for (size_t i = 0; i < num_samples; ++i) {
    class_counts[i].resize(class_values.size(), 0);
  }

  // For each tree loop over OOB samples and count classes
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sample_idx = 0; sample_idx < trees[tree_idx]->getNumSamplesOob(); ++sample_idx) {
      size_t sampleID = trees[tree_idx]->getOobSampleIDs()[sample_idx];
      std::vector<double> counts = trees[tree_idx]->getPredictions()[sample_idx];
      size_t classID = mostFrequentClass(counts, random_number_generator);

      // classID >= counts.size() means all zero
      if (classID < counts.size()) {
        ++class_counts[sampleID][classID];
      }
    }
  }

// Compute majority vote for each sample
  predictions.reserve(num_samples);
  for (size_t i = 0; i < num_samples; ++i) {
    std::vector<double> temp;
    size_t classID = mostFrequentClass(class_counts[i], random_number_generator);

    // classID >= class_counts[i].size() means all zero
    if (classID < class_counts[i].size()) {
      temp.push_back(class_values[classID]);
    } else {
      temp.push_back(NAN);
    }
    predictions.push_back(temp);
  }

  // Compare predictions with true data
  size_t num_missclassifications = 0;
  for (size_t i = 0; i < predictions.size(); ++i) {
    double predicted_value = predictions[i][0];
    if (!std::isnan(predicted_value)) {
      double real_value = data->get(i, dependent_varID);
      if (predicted_value != real_value) {
        ++num_missclassifications;
      }
      ++classification_table[std::make_pair(real_value, predicted_value)];
    }
  }
  overall_prediction_error = (double) num_missclassifications / (double) predictions.size();
}

void ForestProbability::writeOutputInternal() {
  *verbose_out << "Tree type:                         " << "Probability estimation" << std::endl;
}

void ForestProbability::writeConfusionFile() {

  // Open confusion file for writing
  std::string filename = output_prefix + ".confusion";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to confusion file: " + filename + ".");
  }

  // Write confusion to file
  outfile << "Overall OOB prediction error (Fraction missclassified): " << overall_prediction_error << std::endl;
  outfile << std::endl;
  outfile << "Class specific prediction errors:" << std::endl;
  outfile << "           ";
  for (auto& class_value : class_values) {
    outfile << "     " << class_value;
  }
  outfile << std::endl;
  for (auto& predicted_value : class_values) {
    outfile << "predicted " << predicted_value << "     ";
    for (auto& real_value : class_values) {
      size_t value = classification_table[std::make_pair(real_value, predicted_value)];
      outfile << value;
      if (value < 10) {
        outfile << "     ";
      } else if (value < 100) {
        outfile << "    ";
      } else if (value < 1000) {
        outfile << "   ";
      } else if (value < 10000) {
        outfile << "  ";
      } else if (value < 100000) {
        outfile << " ";
      }
    }
    outfile << std::endl;
  }

  outfile.close();
  *verbose_out << "Saved confusion matrix to file " << filename << "." << std::endl;
}

void ForestProbability::writePredictionFile() {

  // Open prediction file for writing
  std::string filename = output_prefix + ".prediction";
  std::ofstream outfile;
  outfile.open(filename, std::ios::out);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to prediction file: " + filename + ".");
  }

  // Write
  outfile << "Class predictions, one sample per row." << std::endl;
  for (auto& class_value : class_values) {
    outfile << class_value << " ";
  }
  outfile << std::endl;
  outfile << std::endl;
  for (size_t i = 0; i < predictions.size(); ++i) {
    for (size_t j = 0; j < predictions[i].size(); ++j) {
      outfile << predictions[i][j] << " ";
    }
    outfile << std::endl;
  }

  *verbose_out << "Saved predictions to file " << filename << "." << std::endl;
}

void ForestProbability::saveToFileInternal(std::ofstream& outfile) {

  // Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

  // Write treetype
  TreeType treetype = TREE_PROBABILITY;
  outfile.write((char*) &treetype, sizeof(treetype));

  // Write class_values
  saveVector1D(class_values, outfile);
}

void ForestProbability::loadFromFileInternal(std::ifstream& infile) {

  // Read number of variables
  size_t num_variables_saved;
  infile.read((char*) &num_variables_saved, sizeof(num_variables_saved));

  // Read treetype
  TreeType treetype;
  infile.read((char*) &treetype, sizeof(treetype));
  if (treetype != TREE_PROBABILITY) {
    throw std::runtime_error("Wrong treetype. Loaded file is not a probability estimation forest.");
  }

  // Read class_values
  readVector1D(class_values, infile);

  for (size_t i = 0; i < num_trees; ++i) {

    // Read data
    std::vector<std::vector<size_t>> child_nodeIDs;
    readVector2D(child_nodeIDs, infile);
    std::vector<size_t> split_varIDs;
    readVector1D(split_varIDs, infile);
    std::vector<double> split_values;
    readVector1D(split_values, infile);

    // Read Terminal node class counts
    std::vector<size_t> terminal_nodes;
    readVector1D(terminal_nodes, infile);
    std::vector<std::vector<double>> terminal_class_counts_vector;
    readVector2D(terminal_class_counts_vector, infile);

    // Convert Terminal node class counts to vector with empty elemtents for non-terminal nodes
    std::vector<std::vector<double>> terminal_class_counts;
    terminal_class_counts.resize(child_nodeIDs.size(), std::vector<double>());
    for (size_t i = 0; i < terminal_nodes.size(); ++i) {
      terminal_class_counts[terminal_nodes[i]] = terminal_class_counts_vector[i];
    }

    // If dependent variable not in test data, change variable IDs accordingly
    if (num_variables_saved > num_variables) {
      for (auto& varID : split_varIDs) {
        if (varID >= dependent_varID) {
          --varID;
        }
      }
    }

    // Create tree
    Tree* tree = new TreeProbability(child_nodeIDs, split_varIDs, split_values, &class_values, &response_classIDs,
        terminal_class_counts);
    trees.push_back(tree);
  }
}

