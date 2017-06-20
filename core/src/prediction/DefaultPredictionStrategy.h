/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_PREDICTIONSTRATEGY_H
#define GRF_PREDICTIONSTRATEGY_H

#include <unordered_map>
#include <vector>

#include "commons/globals.h"
#include "commons/Observations.h"
#include "prediction/Prediction.h"
#include "prediction/PredictionValues.h"

/**
 * A prediction strategy defines how predictions are computed over test samples.
 * This strategy is given a weighted list of training sample IDs that share a leaf
 * with the test sample. To create a more performant strategy, or one that can compute
 * variance estimates, please refer to OptimizedPredictionStrategy.
 */
class DefaultPredictionStrategy {
public:

/**
 * The number of values in a prediction, e.g. 1 for regression
 * or the number of quantiles for quantile forests.
 */
 virtual size_t prediction_length() = 0;

/**
 * Computes a prediction for a single test sample.
 *
 * sample: the ID of the test sample.
 * weights_by_sample: a map from neighboring sample ID, to a weight specifying
 *     how often the sample appeared in the same leaf as the test sample. Note that
 *     these weights are normalized and will sum to 1.
 * observations: the list of observations for all training samples.
 */
virtual std::vector<double> predict(size_t sample,
  const std::unordered_map<size_t, double>& weights_by_sample,
  const Observations& observations) = 0;
};


#endif //GRF_PREDICTIONSTRATEGY_H
