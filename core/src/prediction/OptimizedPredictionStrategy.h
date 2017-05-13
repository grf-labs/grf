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

#ifndef GRADIENTFOREST_OPTIMIZEDPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_OPTIMIZEDPREDICTIONSTRATEGY_H

#include <unordered_map>
#include <vector>

#include "commons/globals.h"
#include "commons/Observations.h"
#include "prediction/Prediction.h"
#include "prediction/PredictionValues.h"

class OptimizedPredictionStrategy {
public:
  virtual size_t prediction_length() = 0;
  virtual Prediction predict(const std::vector<double>& average_prediction_values) = 0;

  virtual Prediction predict_with_variance(size_t sampleID,
                                           const std::vector<std::vector<size_t>>& leaf_sampleIDs,
                                           const Observations& observations,
                                           uint ci_group_size) = 0;

  virtual PredictionValues precompute_prediction_values(
      const std::vector<std::vector<size_t>>& leaf_sampleIDs,
      const Observations& observations) = 0;
};


#endif //GRADIENTFOREST_OPTIMIZEDPREDICTIONSTRATEGY_H
