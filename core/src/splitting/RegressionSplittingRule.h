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

  Authorship: Marvin Wright (wright@imbs.uni-luebeck.de), refactored by
  Julie Tibshirani (jtibs@cs.stanford.edu)
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_REGRESSIONSPLITTINGRULE_H
#define GRADIENTFOREST_REGRESSIONSPLITTINGRULE_H

#include "Tree.h"
#include "SplittingRule.h"
#include <unordered_map>
#include "Data.h"

class RegressionSplittingRule: public SplittingRule {
public:
  RegressionSplittingRule(Data *data);

  ~RegressionSplittingRule();

  bool findBestSplit(size_t nodeID,
                     const std::vector<size_t>& possible_split_varIDs,
                     const std::unordered_map<size_t, double>& labels_by_sampleID,
                     const std::vector<std::vector<size_t>>& sampleIDs,
                     std::vector<size_t>& split_varIDs,
                     std::vector<double>& split_values);

private:
  virtual void findBestSplitValueSmallQ(size_t nodeID,
                                        size_t varID,
                                        double sum_node,
                                        size_t num_samples_node,
                                        double& best_value,
                                        size_t& best_varID,
                                        double& best_decrease,
                                        const std::unordered_map<size_t, double>& responses_by_sampleID,
                                        const std::vector<std::vector<size_t>>& sampleIDs);
  virtual void findBestSplitValueLargeQ(size_t nodeID,
                                        size_t varID,
                                        double sum_node,
                                        size_t num_samples_node,
                                        double& best_value,
                                        size_t& best_varID,
                                        double& best_decrease,
                                        const std::unordered_map<size_t, double>& responses_by_sampleID,
                                        const std::vector<std::vector<size_t>>& sampleIDs);

  Data* data;
  size_t* counter;
  double* sums;

  DISALLOW_COPY_AND_ASSIGN(RegressionSplittingRule);
};


#endif //GRADIENTFOREST_REGRESSIONSPLITTINGRULE_H
