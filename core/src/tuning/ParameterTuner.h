/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#ifndef GRF_PARAMETERTUNER_H
#define GRF_PARAMETERTUNER_H

#include "commons/globals.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestTrainer.h"

class ParameterTuner {

public:
  ParameterTuner(const ForestTrainer& trainer,
                 const ForestPredictor& predictor,
                 uint outcome_index);
  uint tune_min_node_size(Data* data,
                          ForestOptions &options);

private:
  double calculate_mse(const std::vector<Prediction>& predictions,
                       Data* data);

  ForestTrainer trainer;
  ForestPredictor predictor;
  uint outcome_index;
};



#endif //GRF_PARAMETERTUNER_H
