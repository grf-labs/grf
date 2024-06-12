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

#ifndef GRF_MULTIREGRESSIONPREDICTIONSTRATEGY_H
#define GRF_MULTIREGRESSIONPREDICTIONSTRATEGY_H

#include "commons/Data.h"
#include "prediction/OptimizedPredictionStrategy.h"
#include "prediction/PredictionValues.h"
#include "ObjectiveBayesDebiaser.h"

namespace grf {

class MultiRegressionPredictionStrategy final: public OptimizedPredictionStrategy {
public:
  MultiRegressionPredictionStrategy(size_t num_outcomes);

  size_t prediction_value_length() const;

  PredictionValues precompute_prediction_values(const std::vector<std::vector<size_t>>& leaf_samples,
                                                const Data& data) const;

  size_t prediction_length() const;

  std::vector<double> predict(const std::vector<double>& average) const;

  std::vector<double> compute_variance(
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      size_t ci_group_size) const;

  std::vector<std::pair<double, double>> compute_error(
      size_t sample,
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      const Data& data) const;

private:
  size_t num_outcomes;
  size_t num_types;
  size_t weight_index;
};

} // namespace grf

#endif //GRF_MULTIREGRESSIONPREDICTIONSTRATEGY_H
