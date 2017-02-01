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

  Authorship: Julie Tibshirani (jtibs@cs.stanford.edu), loosely based on code
  by Marvin Wright (wright@imbs.uni-luebeck.de)
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_FORESTPREDICTOR_H
#define GRADIENTFOREST_FORESTPREDICTOR_H

#include "RelabelingStrategy.h"
#include "SplittingRule.h"
#include "Prediction.h"
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
  std::map<size_t, std::vector<size_t>> determine_terminal_node_IDs(
      const Forest& forest,
      Data* data,
      bool oob_prediction);

  std::map<size_t, std::vector<size_t>> predictTreesInThread(
      uint thread_idx,
      const Forest& forest,
      Data* prediction_data,
      bool oob_prediction);

  void add_prediction_values(size_t nodeID,
                             const PredictionValues& prediction_values,
                             std::map<std::string, double>& average_prediction_values);

  void add_sample_weights(const std::vector<size_t>& sampleIDs,
                          std::unordered_map<size_t, double>& weights_by_sampleID);

  void normalize_prediction_values(size_t num_leaf_nodes,
                                   std::map<std::string, double>& average_prediction_values);

  void normalize_sample_weights(std::unordered_map<size_t, double>& weights_by_sampleID);

  uint num_threads;
  uint ci_group_size;
  std::vector<uint> thread_ranges;

  std::shared_ptr<PredictionStrategy> prediction_strategy;
};


#endif //GRADIENTFOREST_FORESTPREDICTOR_H
