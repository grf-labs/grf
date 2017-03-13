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

#include "ForestPredictor.h"
#include "OptimizedPredictionCollector.h"
#include "DefaultPredictionCollector.h"
#include "utility.h"

ForestPredictor::ForestPredictor(uint num_threads,
                                 uint ci_group_size,
                                 std::shared_ptr<PredictionStrategy> prediction_strategy) {
  this->prediction_collector = prediction_strategy->requires_leaf_sampleIDs()
      ? std::shared_ptr<PredictionCollector>(new DefaultPredictionCollector(prediction_strategy, ci_group_size))
      : std::shared_ptr<PredictionCollector>(new OptimizedPredictionCollector(prediction_strategy));

  this->num_threads = num_threads == DEFAULT_NUM_THREADS
      ? std::thread::hardware_concurrency()
      : num_threads;
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest, Data* data) {
  std::vector<std::vector<size_t>> leaf_nodes_by_tree = find_leaf_nodes(forest, data, false);
  return prediction_collector->collect_predictions(forest, data, leaf_nodes_by_tree, false);
}

std::vector<Prediction> ForestPredictor::predict_oob(const Forest& forest, Data* data) {
  std::vector<std::vector<size_t>> leaf_nodes_by_tree = find_leaf_nodes(forest, data, true);
  return prediction_collector->collect_predictions(forest, data, leaf_nodes_by_tree, true);
}

std::vector<std::vector<size_t>> ForestPredictor::find_leaf_nodes(
    const Forest& forest,
    Data* data,
    bool oob_prediction) {
  size_t num_trees = forest.get_trees().size();

  std::vector<std::vector<size_t>> leaf_nodes_by_tree;
  leaf_nodes_by_tree.reserve(num_trees);

  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, num_trees - 1, num_threads);

  std::vector<std::future<
      std::vector<std::vector<size_t>>>> futures;
  futures.reserve(num_threads);

  for (uint i = 0; i < num_threads; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_trees_batch = thread_ranges[i + 1] - start_index;
    futures.push_back(std::async(std::launch::async,
                                 &ForestPredictor::find_batch,
                                 this,
                                 start_index,
                                 num_trees_batch,
                                 forest,
                                 data,
                                 oob_prediction));
  }

  for (auto& future : futures) {
    std::vector<std::vector<size_t>> leaf_nodes = future.get();
    leaf_nodes_by_tree.insert(leaf_nodes_by_tree.end(),
                              leaf_nodes.begin(),
                              leaf_nodes.end());
  }

  return leaf_nodes_by_tree;
};

std::vector<std::vector<size_t>> ForestPredictor::find_batch(
    size_t start,
    size_t num_trees,
    const Forest &forest,
    Data *prediction_data,
    bool oob_prediction) {
  std::vector<std::vector<size_t>> all_leaf_node_IDs(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    std::shared_ptr<Tree> tree = forest.get_trees()[start + i];

    const std::vector<size_t> &sampleIDs = oob_prediction ? tree->get_oob_sampleIDs() : std::vector<size_t>();
    std::vector<size_t> leaf_node_IDs = tree->find_leaf_nodeIDs(prediction_data, sampleIDs);
    all_leaf_node_IDs[i] = leaf_node_IDs;
  }

  return all_leaf_node_IDs;
}
