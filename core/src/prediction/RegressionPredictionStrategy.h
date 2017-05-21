/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_REGRESSIONPREDICTIONSTRATEGY_H
#define GRF_REGRESSIONPREDICTIONSTRATEGY_H

#include "commons/Data.h"
#include "commons/Observations.h"
#include "prediction/OptimizedPredictionStrategy.h"
#include "prediction/PredictionValues.h"
#include "ObjectiveBayesDebiaser.h"

class RegressionPredictionStrategy: public OptimizedPredictionStrategy {
public:
  size_t prediction_value_length();

  PredictionValues precompute_prediction_values(const std::vector<std::vector<size_t>>& leaf_sampleIDs,
                                                const Observations& observations);

  size_t prediction_length();

  std::vector<double> predict(const std::vector<double>& average);

  std::vector<double> compute_variance(
      const std::vector<double>& average,
      const PredictionValues& leaf_values,
      uint ci_group_size);

private:
  static const std::size_t OUTCOME;
  ObjectiveBayesDebiaser bayes_debiaser;
};


#endif //GRF_REGRESSIONPREDICTIONSTRATEGY_H
