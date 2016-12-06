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

#include "globals.h"
#include "TreeFactory.h"
#include "Data.h"

class Forest {
public:
  Forest(std::unordered_map<std::string, size_t> observables,
         RelabelingStrategy* relabeling_strategy,
         SplittingRule* splitting_rule);
  virtual ~Forest();

  // Init from c++ main or Rcpp from R
  void initCpp(std::string dependent_variable_name, MemoryMode memory_mode, uint mtry,
      uint num_trees, std::ostream* verbose_out, uint seed, uint num_threads,
      std::string load_forest_filename, uint min_node_size,
      std::string split_select_weights_file, std::vector<std::string>& always_split_variable_names,
      std::string status_variable_name, bool sample_with_replacement,
      bool memory_saving_splitting,
      std::string case_weights_file, double sample_fraction, Data* input_data);
  void init(std::string dependent_variable_name, MemoryMode memory_mode, Data* input_data, uint mtry,
      uint num_trees, uint seed, uint num_threads,
      uint min_node_size, std::string status_variable_name, bool prediction_mode, bool sample_with_replacement,
      bool memory_saving_splitting,
      double sample_fraction);
  virtual void initInternal(std::string status_variable_name) = 0;

  // Grow or predict
  void run(bool verbose);

  // Write results to output files
  void writeOutput();
  virtual void writeOutputInternal() = 0;
  virtual void writeConfusionFile() = 0;
  virtual void writePredictionFile() = 0;

  // Save forest to file
  void saveToFile();
  virtual void saveToFileInternal(std::ofstream& outfile) = 0;

protected:
  void grow();
  void growInternal();

  // Predict using existing tree from file and data as prediction data
  void predict();
  virtual void predictInternal() = 0;

  void computePredictionError();
  virtual void computePredictionErrorInternal() = 0;

  void growTreesInThread(uint thread_idx);
  void predictTreesInThread(uint thread_idx, const Data* prediction_data, bool oob_prediction);

  // Load forest from file
  void loadFromFile(std::string filename);
  virtual void loadFromFileInternal(std::ifstream& infile) = 0;

  // Set split select weights and variables to be always considered for splitting
  void setSplitWeightVector(std::vector<std::vector<double>>& split_select_weights);
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
  bool memory_saving_splitting;
  bool keep_inbag;
  double sample_fraction;

  // Variable to not split at (only dependent_varID for non-survival forests)
  std::vector<size_t> no_split_variables;

  // Multithreading
  uint num_threads;
  std::vector<uint> thread_ranges;
  std::mutex mutex;
  std::condition_variable condition_variable;

  std::vector<TreeFactory*> trees;
  Data* data;
  std::unordered_map<std::string, std::vector<double>> original_observations;

  std::vector<std::vector<double>> predictions;
  double overall_prediction_error;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<size_t> deterministic_varIDs;
  std::vector<size_t> split_select_varIDs;
  std::vector<std::vector<double>> split_select_weights;

  // Bootstrap weights
  std::vector<double> case_weights;

  // Random number generator
  std::mt19937_64 random_number_generator;

  // Computation progress (finished trees)
  size_t progress;

  std::unordered_map<std::string, size_t> observables;
  RelabelingStrategy* relabeling_strategy;
  SplittingRule* splitting_rule;

private:
  DISALLOW_COPY_AND_ASSIGN(Forest);
};

#endif /* FOREST_H_ */
