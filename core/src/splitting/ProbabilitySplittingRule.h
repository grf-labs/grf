/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_PROBABILITYSPLITTINGRULE_H
#define GRADIENTFOREST_PROBABILITYSPLITTINGRULE_H

#include "commons/globals.h"
#include <vector>
#include "commons/Data.h"
#include "splitting/SplittingRule.h"

class ProbabilitySplittingRule: public SplittingRule {
public:
  ProbabilitySplittingRule(Data* data, size_t num_classes);
  ~ProbabilitySplittingRule();

  bool find_best_split(size_t nodeID,
                       const std::vector<size_t>& possible_split_varIDs,
                       const std::unordered_map<size_t, double>& labels_by_sampleID,
                       const std::vector<std::vector<size_t>>& sampleIDs,
                       std::vector<size_t>& split_varIDs,
                       std::vector<double>& split_values);

private:
  void find_best_split_value_small_q(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
                                     size_t num_samples_node,
                                     double& best_value, size_t& best_varID, double& best_decrease,
                                     const std::unordered_map<size_t, double>& labels_by_sampleID,
                                     const std::vector<std::vector<size_t>>& sampleIDs);

  void find_best_split_value_large_q(size_t nodeID, size_t varID, size_t num_classes, size_t* class_counts,
                                     size_t num_samples_node,
                                     double& best_value, size_t& best_varID, double& best_decrease,
                                     const std::unordered_map<size_t, double>& labels_by_sampleID,
                                     const std::vector<std::vector<size_t>>& sampleIDs);

  Data* data;
  size_t num_classes;

  size_t* counter;
  size_t* counter_per_class;

  DISALLOW_COPY_AND_ASSIGN(ProbabilitySplittingRule);
};

#endif //GRADIENTFOREST_PROBABILITYSPLITTINGRULE_H
