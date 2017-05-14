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

#ifndef GRADIENTFOREST_STANDARDPREDICTOR_H
#define GRADIENTFOREST_STANDARDPREDICTOR_H


#include "forest/Forest.h"
#include "prediction/collector/PredictionCollector.h"

class OptimizedPredictionCollector: public PredictionCollector {
public:
  OptimizedPredictionCollector(std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy);

  std::vector<Prediction> collect_predictions(const Forest &forest,
                                              Data *prediction_data,
                                              std::vector<std::vector<size_t>> leaf_nodes_by_tree,
                                              std::vector<std::vector<bool>> trees_by_sample);

private:
  void add_prediction_values(size_t nodeID,
                             const PredictionValues& prediction_values,
                             std::vector<double>& average_prediction_values);

  void normalize_prediction_values(size_t num_leaf_nodes,
                                   std::vector<double>& average_prediction_values);

  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy;
};


#endif //GRADIENTFOREST_STANDARDPREDICTOR_H
