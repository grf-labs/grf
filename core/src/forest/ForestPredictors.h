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

#ifndef GRF_FORESTPREDICTORS_H
#define GRF_FORESTPREDICTORS_H

#include "forest/ForestPredictor.h"

namespace grf {

ForestPredictor instrumental_predictor(uint num_threads);

ForestPredictor multi_causal_predictor(uint num_threads, size_t num_treatments, size_t num_outcomes);

ForestPredictor quantile_predictor(uint num_threads,
                                   const std::vector<double>& quantiles);

ForestPredictor probability_predictor(uint num_threads, size_t num_classes);

ForestPredictor regression_predictor(uint num_threads);

ForestPredictor multi_regression_predictor(uint num_threads, size_t num_outcomes);

ForestPredictor ll_regression_predictor(uint num_threads,
                                        std::vector<double> lambdas,
                                        bool weight_penalty,
                                        std::vector<size_t> linear_correction_variables);

ForestPredictor ll_causal_predictor(uint num_threads,
                                   std::vector<double> lambdas,
                                   bool weight_penalty,
                                   std::vector<size_t> linear_correction_variables);

ForestPredictor survival_predictor(uint num_threads, size_t num_failures, int prediction_type);

ForestPredictor causal_survival_predictor(uint num_threads);

} // namespace grf

#endif //GRF_FORESTPREDICTORS_H
