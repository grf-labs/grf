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

  Author: Julie Tibshirani (jtibs@cs.stanford.edu)
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_PREDICTIONSTRATEGY_H
#define GRADIENTFOREST_PREDICTIONSTRATEGY_H

#include <unordered_map>
#include <vector>

#include "globals.h"
#include "Observations.h"
#include "Prediction.h"
#include "PredictionValues.h"

class PredictionStrategy {
public:
  virtual size_t prediction_length() = 0;
  virtual Prediction predict(const std::map<std::string, double>& average_prediction_values,
                             const std::unordered_map<size_t, double>& weights_by_sampleID,
                             const Observations& observations) = 0;

  virtual Prediction predict_with_variance(const std::vector<std::vector<size_t>>& leaf_sampleIDs,
                                           const Observations& observations,
                                           uint ci_group_size) = 0;

  virtual bool requires_leaf_sampleIDs() = 0;
  virtual PredictionValues precompute_prediction_values(
      const std::vector<std::vector<size_t>>& leaf_sampleIDs,
      const Observations& observations) = 0;
};


#endif //GRADIENTFOREST_PREDICTIONSTRATEGY_H
