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

#ifndef FORESTSURVIVAL_H_
#define FORESTSURVIVAL_H_

#include <iostream>
#include <vector>

#include "globals.h"
#include "Forest.h"
#include "TreeSurvival.h"

class ForestSurvival: public Forest {
public:
  ForestSurvival();
  virtual ~ForestSurvival();

  void loadForest(size_t dependent_varID, size_t num_trees,
      std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      size_t status_varID, std::vector<std::vector<std::vector<double>> >& forest_chf,
      std::vector<double>& unique_timepoints, std::vector<bool>& is_ordered_variable);

  std::vector<std::vector<std::vector<double>>>getChf() {
    std::vector<std::vector<std::vector<double>>> result;
    result.reserve(num_trees);
    for (Tree* tree : trees) {
      TreeSurvival* temp = (TreeSurvival*) tree;
      result.push_back(temp->getChf());
    }
    return result;
  }
  size_t getStatusVarId() const {
    return status_varID;
  }
  const std::vector<double>& getUniqueTimepoints() const {
    return unique_timepoints;
  }

private:
  void initInternal(std::string status_variable_name);
  void growInternal();
  void predictInternal();
  void computePredictionErrorInternal();
  void writeOutputInternal();
  void writeConfusionFile();
  void writePredictionFile();
  void saveToFileInternal(std::ofstream& outfile);
  void loadFromFileInternal(std::ifstream& infile);

  size_t status_varID;
  std::vector<double> unique_timepoints;
  std::vector<size_t> response_timepointIDs;

  DISALLOW_COPY_AND_ASSIGN(ForestSurvival);
};

#endif /* FORESTSURVIVAL_H_ */
