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

#ifndef GRF_FORESTTRAINERS_H
#define GRF_FORESTTRAINERS_H

#include "forest/ForestTrainer.h"

namespace grf {

ForestTrainer instrumental_trainer(double reduced_form_weight,
                                   bool stabilize_splits);

ForestTrainer multi_causal_trainer(size_t num_treatments,
                                   size_t num_outcomes,
                                   bool stabilize_splits,
                                   const std::vector<double>& gradient_weights = {});

ForestTrainer quantile_trainer(const std::vector<double>& quantiles);

ForestTrainer probability_trainer(size_t num_classes);

ForestTrainer regression_trainer();

ForestTrainer multi_regression_trainer(size_t num_outcomes);

ForestTrainer ll_regression_trainer(double split_lambda,
                                   bool weight_penalty,
                                   const std::vector<double>& overall_beta,
                                   size_t ll_split_cutoff,
                                   std::vector<size_t> ll_split_variables);

ForestTrainer survival_trainer();

ForestTrainer causal_survival_trainer(bool stabilize_splits);

} // namespace grf

#endif //GRF_FORESTTRAINERS_H
