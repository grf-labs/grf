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

#ifndef GRF_FORESTPREDICTOR_H
#define GRF_FORESTPREDICTOR_H

#include "relabeling/RelabelingStrategy.h"
#include "splitting/SplittingRule.h"
#include "prediction/Prediction.h"
#include "prediction/collector/TreeTraverser.h"
#include "prediction/collector/PredictionCollector.h"
#include "prediction/collector/SampleWeightComputer.h"
#include "prediction/OptimizedPredictionStrategy.h"
#include "prediction/DefaultPredictionStrategy.h"

#include "tree/Tree.h"
#include "tree/TreeTrainer.h"
#include "forest/Forest.h"

#include <memory>
#include <thread>
#include <future>

class ForestPredictor {
public:
  ForestPredictor(uint num_threads,
                  std::shared_ptr<DefaultPredictionStrategy> strategy);

  ForestPredictor(uint num_threads,
                  uint ci_group_size,
                  std::shared_ptr<OptimizedPredictionStrategy> strategy);

  std::vector<Prediction> predict(const Forest& forest, Data* prediction_data) const;
  std::vector<Prediction> predict_oob(const Forest& forest, Data* original_data) const;

private:
  std::vector<Prediction> predict(const Forest& forest,
                                  Data* data,
                                  bool oob_prediction) const;

private:
  TreeTraverser tree_traverser;
  std::shared_ptr<PredictionCollector> prediction_collector;
};


#endif //GRF_FORESTPREDICTOR_H
