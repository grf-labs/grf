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

#ifndef GRADIENTFOREST_FORESTPREDICTOR_H
#define GRADIENTFOREST_FORESTPREDICTOR_H

#include "RelabelingStrategy.h"
#include "SplittingRule.h"
#include "Prediction.h"
#include "PredictionCollector.h"
#include "PredictionStrategy.h"

#include "Tree.h"
#include "TreeTrainer.h"
#include "Forest.h"

#include <memory>
#include <thread>
#include <future>

class ForestPredictor {
public:
  ForestPredictor(uint num_threads,
                  uint ci_group_size,
                  std::shared_ptr<PredictionStrategy> prediction_strategy);

  std::vector<Prediction> predict(const Forest& forest, Data* prediction_data);
  std::vector<Prediction> predict_oob(const Forest& forest, Data* original_data);

private:
  std::vector<std::vector<size_t>> find_leaf_nodes(
      const Forest &forest,
      Data *data,
      bool oob_prediction);

  std::vector<std::vector<size_t>> find_batch(
      size_t start,
      size_t num_trees,
      const Forest &forest,
      Data *prediction_data,
      bool oob_prediction);


private:
  uint num_threads;
  std::shared_ptr<PredictionCollector> prediction_collector;
};


#endif //GRADIENTFOREST_FORESTPREDICTOR_H
