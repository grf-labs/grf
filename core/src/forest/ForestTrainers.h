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

#ifndef GRADIENTFOREST_FORESTTRAINERS_H
#define GRADIENTFOREST_FORESTTRAINERS_H

#include "forest/ForestTrainer.h"

class ForestTrainers {
public:
  static ForestTrainer instrumental_trainer(Data* data,
                                            size_t outcome_index,
                                            size_t treatment_index,
                                            size_t instrument_index,
                                            double split_regularization);

  static ForestTrainer quantile_trainer(Data* data,
                                        size_t outcome_index,
                                        const std::vector<double>& quantiles);

  static ForestTrainer regression_trainer(Data* data,
                                          size_t outcome_index);

  static ForestTrainer custom_trainer(Data* data,
                                      size_t outcome_index);
};


#endif //GRADIENTFOREST_FORESTTRAINERS_H
