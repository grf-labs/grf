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
  TreeSurvival(std::vector<double>* unique_timepoints, size_t status_varID);

  // Create from loaded forest
  TreeSurvival(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<std::vector<double>> chf, std::vector<double>* unique_timepoints);

  virtual ~TreeSurvival();

  void initInternal();

  void addPrediction(size_t nodeID, size_t sampleID);
  void appendToFileInternal(std::ofstream& file);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);

  const std::vector<std::vector<double> >& getChf() const {
    return chf;
  }

private:
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void createEmptyNodeInternal();

  double computePredictionAccuracyInternal();

  void findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs, size_t num_unique_death_times,
      double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueLogRank(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
      size_t num_unique_death_times, double& best_value, size_t& best_varID, double& best_logrank);
  void findBestSplitValueAUC(size_t nodeID, size_t varID, std::vector<double>& possible_split_values,
      double& best_value, size_t& best_varID, double& best_auc);

  void computeDeathCounts(size_t& num_unique_death_times, size_t nodeID);
  double computeLogRankTest(size_t nodeID, size_t varID, double split_value, size_t num_unique_death_times);
  void computeChildDeathCounts(size_t nodeID, size_t varID, double split_value, size_t* num_samples_left_child);

  // Dirty but fast version for GWAS data
  void findBestSplitValueLogRankGWA(size_t nodeID, size_t varID, double& best_value, size_t& best_varID,
      double& best_logrank);

  double computeAucSplit(size_t nodeID, size_t varID, double split_value);

  void reservePredictionMemory(size_t num_predictions) {
    predictions.resize(num_predictions, std::vector<double>());
    for (auto& sample_vector : predictions) {
      sample_vector.resize(num_timepoints, 0);
    }
  }

  void cleanUpInternal() {
    delete[] num_deaths;
    delete[] num_samples_at_risk;
    delete[] num_deaths_left_child;
    delete[] num_samples_at_risk_left_child;
    delete[] num_deaths_1;
    delete[] num_samples_at_risk_1;
    // num_samples_at_risk_0 and num_deaths_0 are deleted by left_child ones
  }

  size_t status_varID;

  // Unique time points for all individuals (not only this bootstrap), sorted
  std::vector<double>* unique_timepoints;
  size_t num_timepoints;

  // For all terminal nodes CHF for all unique timepoints. For other nodes empty vector.
  std::vector<std::vector<double>> chf;

  // Fields to save to while tree growing
  size_t* num_deaths;
  size_t* num_samples_at_risk;
  size_t* num_deaths_left_child;
  size_t* num_samples_at_risk_left_child;
  size_t* num_deaths_0;
  size_t* num_samples_at_risk_0;
  size_t* num_deaths_1;
  size_t* num_samples_at_risk_1;

  DISALLOW_COPY_AND_ASSIGN(TreeSurvival);
};

#endif /* TREESURVIVAL_H_ */
