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

#include "forest/ForestPredictor.h"
#include "prediction/collector/OptimizedPredictionCollector.h"
#include "prediction/collector/DefaultPredictionCollector.h"
#include "commons/utility.h"

ForestPredictor::ForestPredictor(uint num_threads,
                                 uint ci_group_size,
                                 std::shared_ptr<DefaultPredictionStrategy> strategy) :
    tree_traverser(num_threads) {
  this->prediction_collector = std::shared_ptr<PredictionCollector>(
        new DefaultPredictionCollector(strategy, ci_group_size));
}

ForestPredictor::ForestPredictor(uint num_threads,
                                 uint ci_group_size,
                                 std::shared_ptr<OptimizedPredictionStrategy> strategy) :
    tree_traverser(num_threads) {
  this->prediction_collector = std::shared_ptr<PredictionCollector>(
      new OptimizedPredictionCollector(strategy, ci_group_size));
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest, Data* data) const {
  return predict(forest, data, false);
}

std::vector<Prediction> ForestPredictor::predict_oob(const Forest& forest, Data* data) const {
  return predict(forest, data, true);
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest,
                                                 Data* data,
                                                 bool oob_prediction) const {
  std::vector<std::vector<size_t>> leaf_nodes_by_tree = tree_traverser.get_leaf_nodes(forest, data, oob_prediction);
  std::vector<std::vector<bool>> trees_by_sample = tree_traverser.get_valid_trees_by_sample(forest, data, oob_prediction);

  return prediction_collector->collect_predictions(
      forest, data, leaf_nodes_by_tree, trees_by_sample, oob_prediction);
}
