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

#ifndef GRF_MULTICAUSALPREDICTIONSTRATEGY_H
#define GRF_MULTICAUSALPREDICTIONSTRATEGY_H

#include "commons/Data.h"
#include "prediction/Prediction.h"
#include "prediction/OptimizedPredictionStrategy.h"
#include "prediction/PredictionValues.h"
#include "ObjectiveBayesDebiaser.h"

namespace grf {

/**
 * This prediction strategy is an extension of {@link InstrumentalPredictionStrategy}
 * with a multivariate treatment = instrument, and a vector-valued outcome Y.
 *
 * We estimate a forest weighted (alpha) estimate of beta in
 * Y = const + beta W + e, where Y is a vector-valued outcome and W a k-treatment matrix.
 * The forest weight are optionally multiplied by case sample weights.
 *
 * Some internal details:
 * Estimating an alpha-weighted beta requires storing sufficient statistics
 * for each tree, then averaging these to produce beta. The `predict` method
 * takes these forest-averaged sufficient statistics and produces the final estimate.
 *
 * The internal data structure used for storing these is a std::vector<double>.
 * Five sufficent statistics are needed for the estimate in this strategy,
 * they are detailed in `precompute_prediction_values`.
 *
 * The member variable `num_types` is the total storage count for these five
 * statistics (one double for sample weight, num_treatments *
 * num_treatments for the primed design matrix sum_WW, etc...). They are copied
 * contiguously into a std::vector in the `pre-compute` step, and "re-constructed"
 * with pointers in the `predict` step through Eigen's Map class.
 *
 */
class MultiCausalPredictionStrategy final: public OptimizedPredictionStrategy {
public:

  MultiCausalPredictionStrategy(size_t num_treatments,
                                size_t num_outcomes);

  size_t prediction_value_length() const;
  PredictionValues precompute_prediction_values(
      const std::vector<std::vector<size_t>>& leaf_samples,
      const Data& data) const;

  size_t prediction_length() const;

  std::vector<double> predict(const std::vector<double>& average) const;

  std::vector<double> compute_variance(const std::vector<double>& average,
                                       const PredictionValues& leaf_values,
                                       size_t ci_group_size) const;

  std::vector<std::pair<double, double>> compute_error(
      size_t sample,
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      const Data& data) const;

private:
  size_t num_treatments;
  size_t num_outcomes;
  size_t num_types;
  size_t weight_index;
  size_t Y_index;
  size_t W_index;
  size_t YW_index;
  size_t WW_index;
  ObjectiveBayesDebiaser bayes_debiaser;
};

} // namespace grf

#endif //GRF_MULTICAUSALPREDICTIONSTRATEGY_H
