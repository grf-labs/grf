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

#ifndef TREEREGRESSION_H_
#define TREEREGRESSION_H_

#include "globals.h"
#include "Tree.h"
#include <unordered_map>

class TreeRegression: public Tree {
public:
  TreeRegression();

  // Create from loaded forest
  TreeRegression(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values);

  virtual ~TreeRegression();

  void initInternal();

  virtual bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs,
                     std::unordered_map<size_t, double>& responses_by_sampleID);
  virtual void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                double& best_value, size_t& best_varID, double& best_decrease, std::unordered_map<size_t, double>& responses_by_sampleID);
  virtual void findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node,
                                double& best_value, size_t& best_varID, double& best_decrease, std::unordered_map<size_t, double>& responses_by_sampleID);
  double estimate(size_t nodeID);

  double getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return (split_values[terminal_nodeID]);
  }

  size_t* counter;
  double* sums;

private:
  DISALLOW_COPY_AND_ASSIGN(TreeRegression);
};

#endif /* TREEREGRESSION_H_ */
