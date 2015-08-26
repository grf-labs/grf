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

#ifndef TREESURVIVAL_H_
#define TREESURVIVAL_H_

#include "globals.h"
#include "Tree.h"

class TreeSurvival: public Tree {
public:
  TreeSurvival(std::vector<double>* unique_timepoints, size_t status_varID, std::vector<size_t>* response_timepointIDs);

  // Create from loaded forest
  TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints,
      std::vector<size_t>* response_timepointIDs, std::vector<bool>* is_ordered_variable);

  virtual ~TreeSurvival();

  void initInternal();

  void appendToFileInternal(std::ofstream& file);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);

  const std::vector<std::vector<double> >& getChf() const {
    return chf;
  }

  const std::vector<double>& getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return (chf[terminal_nodeID]);
  }

private:

  void createEmptyNodeInternal();
  double computePredictionAccuracyInternal();

  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);

  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  bool findBestSplitMaxstat(size_t nodeID, std::vector<size_t>& possible_split_varIDs);

  void findBestSplitValueLogRank(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
      double& best_value, size_t& best_varID, double& best_logrank);
  void findBestSplitValueLogRankUnordered(size_t nodeID, size_t varID, std::vector<double>& factor_levels,
      double& best_value, size_t& best_varID, double& best_logrank);
  void findBestSplitValueAUC(size_t nodeID, size_t varID,
      double& best_value, size_t& best_varID, double& best_auc);

  void computeDeathCounts(size_t nodeID);
  void computeChildDeathCounts(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
      size_t* num_samples_right_child, size_t* num_samples_at_risk_right_child, size_t* num_deaths_right_child);

  void computeAucSplit(double time_k, double time_l, double status_k, double status_l, double value_k, double value_l,
      size_t num_splits, std::vector<double>& possible_split_values, double* num_count, double* num_total);

  void findBestSplitValueLogRank(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);
  void findBestSplitValueLogRankUnordered(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);

  void cleanUpInternal() {
    delete[] num_deaths;
    delete[] num_samples_at_risk;
  }

  size_t status_varID;

  // Unique time points for all individuals (not only this bootstrap), sorted
  std::vector<double>* unique_timepoints;
  size_t num_timepoints;
  std::vector<size_t>* response_timepointIDs;

  // For all terminal nodes CHF for all unique timepoints. For other nodes empty vector.
  std::vector<std::vector<double>> chf;

  // Fields to save to while tree growing
  size_t* num_deaths;
  size_t* num_samples_at_risk;

  DISALLOW_COPY_AND_ASSIGN(TreeSurvival);
};

#endif /* TREESURVIVAL_H_ */
