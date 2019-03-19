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

#ifndef GRF_INSTRUMENTALPREDICTIONSTRATEGY_H
#define GRF_INSTRUMENTALPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "commons/DefaultData.h"
#include "prediction/Prediction.h"
#include "prediction/OptimizedPredictionStrategy.h"
#include "prediction/PredictionValues.h"
#include "ObjectiveBayesDebiaser.h"

class InstrumentalPredictionStrategy: public OptimizedPredictionStrategy {
public:
  static const std::size_t OUTCOME;
  static const std::size_t TREATMENT;
  static const std::size_t INSTRUMENT;
  static const std::size_t OUTCOME_INSTRUMENT;
  static const std::size_t TREATMENT_INSTRUMENT;

  size_t prediction_value_length();
  PredictionValues precompute_prediction_values(
      const std::vector<std::vector<size_t>>& leaf_samples,
      const Data* data);

  size_t prediction_length();

  std::vector<double> predict(const std::vector<double>& average);

  std::vector<double> compute_variance(const std::vector<double>& average,
                          const PredictionValues& leaf_values,
                          size_t ci_group_size);

  std::vector<std::pair<double, double>>  compute_error(
      size_t sample,
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      const Data* data);

private:
  ObjectiveBayesDebiaser bayes_debiaser;
};

#endif //GRF_INSTRUMENTALPREDICTIONSTRATEGY_H
