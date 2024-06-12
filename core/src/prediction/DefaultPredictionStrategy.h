/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#ifndef GRF_DEFAULTPREDICTIONSTRATEGY_H
#define GRF_DEFAULTPREDICTIONSTRATEGY_H

#include <unordered_map>
#include <vector>

#include "commons/globals.h"
#include "commons/Data.h"
#include "prediction/Prediction.h"
#include "prediction/PredictionValues.h"

namespace grf {

/**
 * A prediction strategy defines how predictions are computed over test samples.
 * This strategy is given a weighted list of training sample IDs that share a leaf
 * with the test sample. To create a more performant strategy, or one that can compute
 * variance estimates, please refer to OptimizedPredictionStrategy.
 */
class DefaultPredictionStrategy {
public:

  virtual ~DefaultPredictionStrategy() = default;

  /**
   * The number of values in a prediction, e.g. 1 for regression
   * or the number of quantiles for quantile forests.
   */
  virtual size_t prediction_length() const = 0;

  /**
   * Computes a prediction for a single test sample.
   *
   * sample: the ID of the test sample.
   * weights_by_sample: a map from neighboring sample ID, to a weight specifying
   *     how often the sample appeared in the same leaf as the test sample. Note that
   *     these weights are normalized and will sum to 1.
   * train_data: the training data matrix.
   * data: the test data matrix. Note that in the case of OOB prediction, this could
   *     be the same as the training matrix.
   */
  virtual std::vector<double> predict(size_t sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Data& train_data,
    const Data& data) const = 0;

  /**
   * Computes a prediction variance estimate for a single test sample.
   *
   * sample: the ID of the test sample.
   * samples_by_tree: vector of samples in the same leaf as the test point,
   *    for each tree
   * weights_by_sampleID: a map from neighboring sample ID, to a weight specifying
   *     how often the sample appeared in the same leaf as the test sample. Note that
   *     these weights are normalized and will sum to 1.
   * train_data: the training data matrix.
   * data: the test data matrix. Note that in the case of OOB prediction, this could
   *     be the same as the training matrix.
   * ci_group_size: the size of the tree groups used to train the forest. This
   *     parameter is used when computing within vs. across group variance.
   */
  virtual std::vector<double> compute_variance(
      size_t sample,
      const std::vector<std::vector<size_t>>& samples_by_tree,
      const std::unordered_map<size_t, double>& weights_by_sampleID,
      const Data& train_data,
      const Data& data,
      size_t ci_group_size) const = 0;
};

} // namespace grf

#endif //GRF_DEFAULTPREDICTIONSTRATEGY_H
