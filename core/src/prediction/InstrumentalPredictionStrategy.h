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

#ifndef GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H
#define GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "commons/Data.h"
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

  size_t prediction_length();
  std::vector<double> predict(const std::vector<double>& average);

  std::vector<double> compute_variance(const std::vector<double>& average,
                          const std::vector<std::vector<double>>& leaf_values,
                          uint ci_group_size);

  PredictionValues precompute_prediction_values(
      const std::vector<std::vector<size_t>>& leaf_sampleIDs,
      const Observations& observations);

private:
  ObjectiveBayesDebiaser bayes_debiaser;
};

#endif //GRADIENTFOREST_INSTRUMENTALPREDICTIONSTRATEGY_H
