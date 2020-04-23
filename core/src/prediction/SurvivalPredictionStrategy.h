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
#include "prediction/OptimizedPredictionStrategy.h"
#include "prediction/PredictionValues.h"

namespace grf {

class SurvivalPredictionStrategy final: public OptimizedPredictionStrategy {
public:

  /**
   * Compute Kaplan-Meier estimates of the survival function.
   *
   * (Variance estimates and other prediction statistics are currently not implemented).
   *
   * num_failures: the count of failures in the training data.
   */
  SurvivalPredictionStrategy(size_t num_failures);

  size_t prediction_length() const;

  std::vector<double> predict(const std::vector<double>& average) const;

  std::vector<double> compute_variance(
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      size_t ci_group_size) const;

  size_t prediction_value_length() const;

  /**
   * Compute Kaplan-Meier estimates of the survival function S(t).
   *
   * The prediction for each sample will be a vector of length `num_failures` with each entry
   * corresponding to an element of S(i) for that individual.
   */
  PredictionValues precompute_prediction_values(const std::vector<std::vector<size_t>>& leaf_samples,
                                                const Data& data) const;

  std::vector<std::pair<double, double>> compute_error(
      size_t sample,
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      const Data& data) const;

private:
  size_t num_failures;
};

} // namespace grf

#endif //GRF_SURVIVALPREDICTIONSTRATEGY_H
