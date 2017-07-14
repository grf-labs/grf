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
                                 std::shared_ptr<DefaultPredictionStrategy> strategy) {
  this->prediction_collector = std::shared_ptr<PredictionCollector>(
        new DefaultPredictionCollector(strategy));
  this->num_threads = num_threads == DEFAULT_NUM_THREADS
      ? std::thread::hardware_concurrency()
      : num_threads;
}

ForestPredictor::ForestPredictor(uint num_threads,
                                 uint ci_group_size,
                                 std::shared_ptr<OptimizedPredictionStrategy> strategy) {
  this->prediction_collector = std::shared_ptr<PredictionCollector>(
      new OptimizedPredictionCollector(strategy, ci_group_size));
  this->num_threads = num_threads == DEFAULT_NUM_THREADS
                      ? std::thread::hardware_concurrency()
                      : num_threads;
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest, Data* data) {
  return predict(forest, data, false);
}

std::vector<Prediction> ForestPredictor::predict_oob(const Forest& forest, Data* data) {
  return predict(forest, data, true);
}

std::vector<Prediction> ForestPredictor::predict(const Forest& forest,
                                                 Data* data,
                                                 bool oob_prediction) {
  std::vector<std::vector<size_t>> leaf_nodes_by_tree = find_leaf_nodes(forest, data, oob_prediction);
  std::vector<std::vector<bool>> trees_by_sample = oob_prediction
          ? get_trees_by_sample(forest, data)
          : std::vector<std::vector<bool>>();
  return prediction_collector->collect_predictions(forest, data, leaf_nodes_by_tree, trees_by_sample);
}

std::vector<std::vector<bool>> ForestPredictor::get_trees_by_sample(const Forest &forest,
                                                                    Data *data) {
  size_t num_trees = forest.get_trees().size();
  size_t num_samples = data->get_num_rows();

  std::vector<std::vector<bool>> result(num_samples, std::vector<bool>(num_trees));

  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sample : forest.get_trees()[tree_idx]->get_oob_samples()) {
      result[sample][tree_idx] = true;
    }
  }
  return result;
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
  futures.reserve(thread_ranges.size());

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
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
  std::vector<std::vector<size_t>> all_leaf_nodes(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    std::shared_ptr<Tree> tree = forest.get_trees()[start + i];

    const std::vector<size_t>& samples = oob_prediction ? tree->get_oob_samples() : std::vector<size_t>();
    std::vector<size_t> leaf_nodes = tree->find_leaf_nodes(prediction_data, samples);
    all_leaf_nodes[i] = leaf_nodes;
  }

  return all_leaf_nodes;
}
