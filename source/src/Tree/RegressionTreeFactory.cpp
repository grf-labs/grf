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
#include <iostream>
#include <vector>
#include <set>

#include "utility.h"
#include "RegressionTreeFactory.h"

RegressionTreeFactory::RegressionTreeFactory() {}

RegressionTreeFactory::RegressionTreeFactory(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values) :
    TreeFactory(child_nodeIDs, split_varIDs, split_values) {
}

double RegressionTreeFactory::estimate(size_t nodeID) {
  double sum_responses_in_node = 0;
  size_t num_samples_in_node = sampleIDs[nodeID].size();
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    sum_responses_in_node += data->get(sampleIDs[nodeID][i], dependent_varID);
  }
  return (sum_responses_in_node / (double) num_samples_in_node);
}

bool RegressionTreeFactory::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

// Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  // Check if node is pure and set split_value to estimate and stop if pure
  bool pure = true;
  double pure_value = 0;
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    double value = data->get(sampleIDs[nodeID][i], dependent_varID);
    if (i != 0 && value != pure_value) {
      pure = false;
      break;
    }
    pure_value = value;
  }
  if (pure) {
    split_values[nodeID] = pure_value;
    return true;
  }

  std::unordered_map<size_t, double> responses_by_sampleID;
  for (size_t sampleID: sampleIDs[nodeID]) {
    responses_by_sampleID[sampleID] = data->get(sampleID, dependent_varID);
  }

  // Find best split, stop if no decrease of impurity
  RegressionSplittingRule *splittingRule = new RegressionSplittingRule(data);
  bool stop = splittingRule->findBestSplit(nodeID,
                                           possible_split_varIDs,
                                           responses_by_sampleID,
                                           sampleIDs,
                                           split_varIDs,
                                           split_values);

  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}