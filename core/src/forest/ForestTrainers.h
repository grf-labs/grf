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

namespace grf {

ForestTrainer instrumental_trainer(double reduced_form_weight,
                                   bool stabilize_splits);

ForestTrainer quantile_trainer(const std::vector<double>& quantiles);

ForestTrainer regression_trainer();

ForestTrainer custom_trainer();

} // namespace grf

#endif //GRF_FORESTTRAINERS_H
