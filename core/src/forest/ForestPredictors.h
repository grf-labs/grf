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

#ifndef GRF_FORESTPREDICTORS_H
#define GRF_FORESTPREDICTORS_H

#include "forest/ForestPredictor.h"

class ForestPredictors {
public:
  static ForestPredictor custom_predictor(uint num_threads);

  static ForestPredictor instrumental_predictor(uint num_threads,
                                                uint ci_group_size);

  static ForestPredictor quantile_predictor(uint num_threads,
                                            const std::vector<double>& quantiles);

  static ForestPredictor regression_predictor(uint num_threads,
                                              uint ci_group_size);

  static ForestPredictor local_linear_predictor(uint num_threads,
                                                uint ci_group_size,
                                                const Data* original_data,
                                                const Data* test_data,
                                                std::vector<double> lambdas,
                                                bool weighted_penalty,
                                                std::vector<size_t> linear_correction_variables);

private:
  static uint get_num_threads(uint provided_num_threads);
};


#endif //GRF_FORESTPREDICTORS_H
