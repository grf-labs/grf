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

#include <stdexcept>

#include "forest/ForestPredictor.h"
#include "prediction/collector/OptimizedPredictionCollector.h"
#include "prediction/collector/DefaultPredictionCollector.h"
#include "commons/utility.h"

namespace grf {

ForestPredictor::ForestPredictor(uint num_threads,
                                 std::unique_ptr<DefaultPredictionStrategy> strategy) :
    tree_traverser(num_threads) {
  this->prediction_collector = std::unique_ptr<PredictionCollector>(
        new DefaultPredictionCollector(std::move(strategy), num_threads));
}

ForestPredictor::ForestPredictor(uint num_threads,
                                 std::unique_ptr<OptimizedPredictionStrategy> strategy) :
    tree_traverser(num_threads) {
  this->prediction_collector = std::unique_ptr<PredictionCollector>(
      new OptimizedPredictionCollector(std::move(strategy), num_threads));
}


std::vector<Prediction> ForestPredictor::predict(const Forest& forest,
                                                 const Data& train_data,
                                                 const Data& data,
                                                 bool estimate_variance) const {
  return predict(forest, train_data, data, estimate_variance, false);
}

std::vector<Prediction> ForestPredictor::predict_oob(const Forest& forest,
                                                     const Data& data,
                                                     bool estimate_variance) const {
  return predict(forest, data, data, estimate_variance, true);
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest,
                                                 const Data& train_data,
                                                 const Data& data,
                                                 bool estimate_variance,
                                                 bool oob_prediction) const {
    // debug information
    std::cout << "Entering ForestPredictor::predict" << std::endl;
    std::cout << "Number of trees: " << forest.get_trees().size() << std::endl;
    std::cout << "Train data dimensions: " << train_data.get_num_rows() << "x" << train_data.get_num_cols() << std::endl;
    std::cout << "Test data dimensions: " << data.get_num_rows() << "x" << data.get_num_cols() << std::endl;

    if (estimate_variance && forest.get_ci_group_size() <= 1) {
        throw std::runtime_error("To estimate variance during prediction, the forest must"
           " be trained with ci_group_size greater than 1.");
    }

    std::vector<std::vector<size_t>> leaf_nodes_by_tree = tree_traverser.get_leaf_nodes(forest, data, oob_prediction);
    // debug
    std::cout << "Leaf nodes by tree size: " << leaf_nodes_by_tree.size() << std::endl;

    std::vector<std::vector<bool>> trees_by_sample = tree_traverser.get_valid_trees_by_sample(forest, data, oob_prediction);
    // debug
    std::cout << "Trees by sample size: " << trees_by_sample.size() << std::endl;

    return prediction_collector->collect_predictions(forest, train_data, data,
        leaf_nodes_by_tree, trees_by_sample,
        estimate_variance, oob_prediction);
}

} // namespace grf
