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

#ifndef GRF_CUSTOMPREDICTIONSTRATEGY_H
#define GRF_CUSTOMPREDICTIONSTRATEGY_H


#include "Eigen/Dense"
#include "DefaultPredictionStrategy.h"

class CustomPredictionStrategy: public DefaultPredictionStrategy {
public:
  // Add more observables here as needed.
  static const std::size_t OUTCOME;

  size_t prediction_value_length();

  size_t prediction_length();
  std::vector<double> predict(size_t sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Observations& observations);

  std::vector<double> compute_variance(
          std::vector<std::vector<size_t>> samples_by_tree,
          uint ci_group_size,
          size_t sampleID,
          std::unordered_map<size_t, double> weights_by_sampleID,
          const Observations& observations,
          const PredictionValues& leaf_values);

  std::vector<double> compute_debiased_error(
          size_t sample,
          const PredictionValues& leaf_values,
          const Observations& observations);

};


#endif //GRF_CUSTOMPREDICTIONSTRATEGY_H
