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

#ifndef GRF_LLCAUSALPREDICTIONSTRATEGY_H
#define GRF_LLCAUSALPREDICTIONSTRATEGY_H


#include <cstddef>
#include <unordered_map>
#include "Eigen/Dense"
#include "commons/Data.h"
#include "prediction/Prediction.h"
#include "prediction/DefaultPredictionStrategy.h"
#include "prediction/PredictionValues.h"
#include "ObjectiveBayesDebiaser.h"

namespace grf {

class LLCausalPredictionStrategy final: public DefaultPredictionStrategy {
public:
    LLCausalPredictionStrategy(std::vector<double> lambdas,
                               bool weight_penalty,
                               std::vector<size_t> linear_correction_variables);

    size_t prediction_length() const;

    std::vector<double> predict(size_t sampleID,
                                const std::unordered_map<size_t, double>& weights_by_sampleID,
                                const Data& original_data,
                                const Data& test_data) const;

    std::vector<double> compute_variance(
            size_t sampleID,
            const std::vector<std::vector<size_t>>& samples_by_tree,
            const std::unordered_map<size_t, double>& weights_by_sampleID,
            const Data& train_data,
            const Data& data,
            size_t ci_group_size) const;

private:
    std::vector<double> lambdas;
    bool weight_penalty;
    std::vector<size_t> linear_correction_variables;
    ObjectiveBayesDebiaser bayes_debiaser;
};

} // namespace grf

#endif //GRF_LLCAUSALPREDICTIONSTRATEGY_H
