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

#ifndef FOREST_H_
#define FOREST_H_

#include <vector>
#include <iostream>
#include <random>
#include <mutex>
#include <condition_variable>

#include "globals.h"
#include "Tree.h"
#include "Data.h"

class Forest {
public:
  Forest();
  virtual ~Forest();

  // Init from c++ main or Rcpp from R
  void initCpp(std::string dependent_variable_name, MemoryMode memory_mode, std::string input_file, uint mtry,
      std::string output_prefix, uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
      std::string load_forest_filename, ImportanceMode importance_mode, uint min_node_size,
      std::string split_select_weights_file, std::vector<std::string>& always_split_variable_names,
      std::string status_variable_name, bool sample_with_replacement);
  void initR(std::string dependent_variable_name, MemoryMode memory_mode, Data* input_data, uint mtry, uint num_trees,
      std::ostream* verbose_out, uint seed, uint num_threads, ImportanceMode importance_mode, uint min_node_size,
      std::vector<double>& split_select_weights, std::vector<std::string>& always_split_variable_names,
      std::string status_variable_name, bool prediction_mode, bool sample_with_replacement);
  void init(std::string dependent_variable_name, MemoryMode memory_mode, Data* input_data, uint mtry,
      std::string output_prefix, uint num_trees, uint seed, uint num_threads, ImportanceMode importance_mode,
      uint min_node_size, std::string status_variable_name, bool prediction_mode, bool sample_with_replacement);
  virtual void initInternal(std::string status_variable_name) = 0;

  // Grow or predict
  void run(bool verbose);

  // Write results to output files
  void writeOutput();
  virtual void writeOutputInternal() = 0;
  virtual void writeConfusionFile() = 0;
  virtual void writePredictionFile() = 0;
  void writeImportanceFile();

  // Save forest to file
  void saveToFile();
  virtual void saveToFileInternal(std::ofstream& outfile) = 0;

  std::vector<std::vector<std::vector<size_t>>>getChildNodeIDs() {
    std::vector<std::vector<std::vector<size_t>>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getChildNodeIDs());
    }
    return result;
  }
  std::vector<std::vector<size_t>> getSplitVarIDs() {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getSplitVarIDs());
    }
    return result;
  }
  std::vector<std::vector<double>> getSplitValues() {
    std::vector<std::vector<double>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getSplitValues());
    }
    return result;
  }
  const std::vector<double>& getVariableImportance() const {
    return variable_importance;
  }
  double getOverallPredictionError() const {
    return overall_prediction_error;
  }
  const std::vector<std::vector<double> >& getPredictions() const {
    return predictions;
  }
  size_t getDependentVarId() const {
    return dependent_varID;
  }
  size_t getNumTrees() const {
    return num_trees;
  }
  uint getMtry() const
  {
    return mtry;
  }
  uint getMinNodeSize() const
  {
    return min_node_size;
  }
  size_t getNumIndependentVariables() const
  {
    return num_independent_variables;
  }

protected:
  void grow();
  virtual void growInternal() = 0;

  // Predict using existing tree from file and data as prediction data
  void predict();
  virtual void predictInternal() = 0;

  void computePredictionError();
  virtual void computePredictionErrorInternal() = 0;

  void computeGiniImportance();
  void computePermutationImportance();

  // Multithreading methods for growing/prediction/importance, called by each thread
  void growTreesInThread(uint thread_idx);
  void predictTreesInThread(uint thread_idx, const Data* prediction_data, bool oob_prediction);
  void computeTreePermutationImportanceInThread(uint thread_idx);

  // Load forest from file
  void loadFromFile(std::string filename);
  virtual void loadFromFileInternal(std::ifstream& infile) = 0;

  // Set split select weights and variables to be always considered for splitting
  void setSplitWeightVector(std::vector<double>& split_select_weights);
  void setAlwaysSplitVariables(std::vector<std::string>& always_split_variable_names);

  // Show progress every few seconds
  void showProgress(std::string operation);

  // Verbose output stream, cout if verbose==true, logfile if not
  std::ostream* verbose_out;

  size_t num_trees;
  uint mtry;
  uint min_node_size;
  size_t num_variables;
  size_t num_independent_variables;
  uint seed;
  size_t dependent_varID;
  size_t num_samples;
  bool prediction_mode;
  MemoryMode memory_mode;
  bool sample_with_replacement;

  // Variable to not split at (only dependent_varID for non-survival forests)
  std::vector<size_t> no_split_variables;

  // Multithreading
  uint num_threads;
  std::vector<uint> thread_ranges;
  std::mutex mutex;
  std::condition_variable condition_variable;

  std::vector<Tree*> trees;
  Data* data;

  std::vector<std::vector<double>> predictions;
  double overall_prediction_error;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<size_t> deterministic_varIDs;
  std::vector<size_t> split_select_varIDs;
  std::vector<double> split_select_weights;

  // Random number generator
  std::mt19937_64 random_number_generator;

  std::string output_prefix;
  ImportanceMode importance_mode;

  // Variable importance for all variables in forest
  std::vector<double> variable_importance;

  // Computation progress (finished trees)
  size_t progress;

private:
  DISALLOW_COPY_AND_ASSIGN(Forest);
};

#endif /* FOREST_H_ */
