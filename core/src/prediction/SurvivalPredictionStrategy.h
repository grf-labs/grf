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

#ifndef GRF_SURVIVALPREDICTIONSTRATEGY_H
#define GRF_SURVIVALPREDICTIONSTRATEGY_H

#include "commons/DefaultData.h"
#include "commons/Data.h"
#include "prediction/DefaultPredictionStrategy.h"
#include "prediction/PredictionValues.h"

namespace grf {

class SurvivalPredictionStrategy final: public DefaultPredictionStrategy {
public:

  /**
   * Compute Kaplan-Meier estimates of the survival function.
   *
   * This estimate is weighted by the random forest weights (alpha) and
   * the sample weights.
   *
   * (Variance and error estimates are not supported).
   *
   * num_failures: the count of failures in the training data.
   * The event times retrieved from data.get_outcome(sample) will always be
   * integers in the range 0, ..., num_failures.
   *
   */
  SurvivalPredictionStrategy(size_t num_failures);

  size_t prediction_length() const;

  std::vector<double> predict(size_t prediction_sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Data& train_data,
    const Data& data) const;

  std::vector<double> compute_variance(
    size_t sample,
    const std::vector<std::vector<size_t>>& samples_by_tree,
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    const Data& train_data,
    const Data& data,
    size_t ci_group_size) const;

private:
  size_t num_failures;
};

} // namespace grf

#endif //GRF_SURVIVALPREDICTIONSTRATEGY_H
