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

#ifndef GRF_LLCAUSALPREDICTIONSTRATEGY_H
#define GRF_LLCAUSALPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "commons/DefaultData.h"
#include "prediction/Prediction.h"
#include "prediction/DefaultPredictionStrategy.h"
#include "prediction/PredictionValues.h"

class LLCausalPredictionStrategy: public DefaultPredictionStrategy {
public:
  LLCausalPredictionStrategy(std::vector<double> lambdas,
                           bool use_unweighted_penalty,
                           std::vector<size_t> linear_correction_variables);

  static const std::size_t OUTCOME;
  static const std::size_t TREATMENT;

  size_t prediction_value_length();

  size_t prediction_length();
  std::vector<double> predict(size_t sampleID,
                              const std::unordered_map<size_t, double>& weights_by_sampleID,
                              const Data *original_data,
                              const Data *test_data);

  std::vector<double> compute_variance(
          size_t sampleID,
          std::vector<std::vector<size_t>> samples_by_tree,
          std::unordered_map<size_t, double> weights_by_sampleID,
          const Data* train_data,
          const Data* data,
          size_t ci_group_size);

private:
  const Data *original_data;
  const Data *test_data;
  std::vector<double> lambdas;
  bool use_unweighted_penalty;
  std::vector<size_t> linear_correction_variables;
};

#endif //GRF_LLCAUSALPREDICTIONSTRATEGY_H
