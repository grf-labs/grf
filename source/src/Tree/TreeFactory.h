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

#ifndef TREE_H_
#define TREE_H_

#include <vector>
#include <random>
#include <iostream>

#include "globals.h"
#include "Data.h"
#include "SplittingRule.h"
#include "RelabelingStrategy.h"

class TreeFactory {
public:
  TreeFactory(RelabelingStrategy* relabeling_strategy,
              SplittingRule* splitting_rule);

  // Create from loaded forest
  TreeFactory(std::vector<std::vector<size_t>>& child_nodeIDs,
              std::vector<size_t>& split_varIDs,
              std::vector<double>& split_values,
              std::vector<std::vector<size_t>> sampleIDs,
              RelabelingStrategy* relabeling_strategy,
              SplittingRule* splitting_rule);

  virtual ~TreeFactory();

  void init(Data* data, uint mtry, size_t dependent_varID, size_t num_samples, uint seed,
      std::vector<size_t>* deterministic_varIDs, std::vector<size_t>* split_select_varIDs,
      std::vector<double>* split_select_weights, uint min_node_size,
      std::vector<size_t>* no_split_variables, bool sample_with_replacement,
      std::vector<double>* case_weights, bool keep_inbag,
      double sample_fraction);

  void grow();

  void predict(const Data* prediction_data, bool oob_prediction);

  void appendToFile(std::ofstream& file);

  const std::vector<std::vector<size_t> >& getChildNodeIDs() const {
    return child_nodeIDs;
  }

  const std::vector<std::vector<size_t>>& get_sampleIDs() const {
    return sampleIDs;
  }

  const std::vector<double>& getSplitValues() const {
    return split_values;
  }

  const std::vector<size_t>& getSplitVarIDs() const {
    return split_varIDs;
  }

  const std::vector<size_t>& getOobSampleIDs() const {
    return oob_sampleIDs;
  }
  size_t getNumSamplesOob() const {
    return num_samples_oob;
  }

  const std::vector<size_t>& getInbagCounts() const {
    return inbag_counts;
  }

  std::vector<size_t> get_neighboring_samples(size_t sampleID);

protected:
  void createPossibleSplitVarSubset(std::vector<size_t>& result);

  bool splitNode(size_t nodeID);
  bool splitNodeInternal(size_t nodeID, std::vector<size_t> &possible_split_varIDs);

  void createEmptyNode();

  void bootstrap();
  void bootstrapWithoutReplacement();

  void bootstrapWeighted();
  void bootstrapWithoutReplacementWeighted();

  RelabelingStrategy* relabeling_strategy;
  SplittingRule* splitting_rule;

  Data* data;
  size_t dependent_varID;
  size_t num_samples;

  uint min_node_size;

  std::vector<size_t> split_varIDs;
  std::vector<double> split_values;
  std::vector<std::vector<size_t>> child_nodeIDs;
  std::vector<std::vector<size_t>> sampleIDs;

  // When growing here the OOB set is used
  // Terminal nodeIDs for prediction samples
  std::vector<size_t> prediction_terminal_nodeIDs;

  // Variables related to bootstrapping
  bool sample_with_replacement;
  double sample_fraction;
  bool keep_inbag;
  std::vector<size_t> inbag_counts;
  std::vector<double>* case_weights;
  std::vector<size_t> oob_sampleIDs;
  size_t num_samples_oob;
  std::mt19937_64 random_number_generator;

  // Variables related to variable selection
  std::vector<size_t>* deterministic_varIDs;
  std::vector<size_t>* split_select_varIDs;
  std::vector<double>* split_select_weights;
  std::vector<size_t>* no_split_variables;

  uint mtry;

private:
  DISALLOW_COPY_AND_ASSIGN(TreeFactory);
};

#endif /* TREE_H_ */
