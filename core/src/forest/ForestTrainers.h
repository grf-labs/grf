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

#ifndef GRF_FORESTTRAINERS_H
#define GRF_FORESTTRAINERS_H

#include "forest/ForestTrainer.h"

class ForestTrainers {
public:
  static ForestTrainer instrumental_trainer(size_t outcome_index,
                                            size_t treatment_index,
                                            size_t instrument_index,
                                            double split_regularization,
                                            double alpha,
                                            const ForestOptions& options);

  static ForestTrainer quantile_trainer(size_t outcome_index,
                                        const std::vector<double>& quantiles,
                                        double alpha,
                                        const ForestOptions& options);

  static ForestTrainer regression_trainer(size_t outcome_index,
                                          double alpha,
                                          const ForestOptions& options);

  static ForestTrainer custom_trainer(size_t outcome_index,
                                      double alpha,
                                      const ForestOptions& options);

  static ForestTrainer regularized_regression_trainer(size_t outcome_index,
                                                      double lambda,
                                                      bool downweight_penalty,
                                                      const ForestOptions& options);

  static ForestTrainer regularized_instrumental_trainer(size_t outcome_index,
                                                        size_t treatment_index,
                                                        size_t instrument_index,
                                                        double split_regularization,
                                                        double lambda,
                                                        bool downweight_penalty,
                                                        const ForestOptions& options);
};


#endif //GRF_FORESTTRAINERS_H
