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

#include <future>
#include <stdexcept>

#include "prediction/collector/DefaultPredictionCollector.h"
#include "commons/utility.h"

namespace grf {

DefaultPredictionCollector::DefaultPredictionCollector(std::unique_ptr<DefaultPredictionStrategy> strategy,
                                                       uint num_threads):
    strategy(std::move(strategy)), num_threads(num_threads) {}

std::vector<Prediction> DefaultPredictionCollector::collect_predictions(
    const Forest& forest,
    const Data& train_data,
    const Data& data,
    const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
    const std::vector<std::vector<bool>>& valid_trees_by_sample,
    bool estimate_variance,
    bool estimate_error) const {

  size_t num_samples = data.get_num_rows();
  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, num_samples - 1, num_threads);

  std::vector<std::future<std::vector<Prediction>>> futures;
  futures.reserve(thread_ranges.size());

  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_samples_batch = thread_ranges[i + 1] - start_index;

    futures.push_back(std::async(std::launch::async,
                                 &DefaultPredictionCollector::collect_predictions_batch,
                                 this,
                                 std::ref(forest),
                                 std::ref(train_data),
                                 std::ref(data),
                                 std::ref(leaf_nodes_by_tree),
                                 std::ref(valid_trees_by_sample),
                                 estimate_variance,
                                 start_index,
                                 num_samples_batch));
  }

  for (auto& future : futures) {
    std::vector<Prediction> thread_predictions = future.get();
    predictions.insert(predictions.end(),
                       std::make_move_iterator(thread_predictions.begin()),
                       std::make_move_iterator(thread_predictions.end()));
  }

  return predictions;
}

std::vector<Prediction> DefaultPredictionCollector::collect_predictions_batch(
    const Forest& forest,
    const Data& train_data,
    const Data& data,
    const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
    const std::vector<std::vector<bool>>& valid_trees_by_sample,
    bool estimate_variance,
    size_t start,
    size_t num_samples) const {
  size_t num_trees = forest.get_trees().size();
  bool record_leaf_samples = estimate_variance;

  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  for (size_t sample = start; sample < num_samples + start; ++sample) {
    std::unordered_map<size_t, double> weights_by_sample = weight_computer.compute_weights(
        sample, forest, leaf_nodes_by_tree, valid_trees_by_sample);
    std::vector<std::vector<size_t>> samples_by_tree;

    // If this sample has no neighbors, then return placeholder predictions. Note
    // that this can only occur when honesty is enabled, and is expected to be rare.
    if (weights_by_sample.empty()) {
      std::vector<double> nan(strategy->prediction_length(), NAN);
      predictions.emplace_back(nan, nan, nan, nan);
      continue;
    }

    if (record_leaf_samples) {
      samples_by_tree.resize(num_trees);

      for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
        if (!valid_trees_by_sample[sample][tree_index]) {
          continue;
        }
        const std::vector<size_t>& leaf_nodes = leaf_nodes_by_tree.at(tree_index);
        size_t node = leaf_nodes.at(sample);

        const std::unique_ptr<Tree>& tree = forest.get_trees()[tree_index];
        const std::vector<std::vector<size_t>>& leaf_samples = tree->get_leaf_samples();
        samples_by_tree.push_back(leaf_samples.at(node));
      }
    }

    std::vector<double> point_prediction = strategy->predict(sample, weights_by_sample, train_data, data);
    std::vector<double> variance = estimate_variance
        ? strategy->compute_variance(sample, samples_by_tree, weights_by_sample, train_data, data, forest.get_ci_group_size())
        : std::vector<double>();

    Prediction prediction(point_prediction, variance, {}, {});
    validate_prediction(sample, point_prediction);
    predictions.push_back(prediction);
  }

  return predictions;
}

void DefaultPredictionCollector::validate_prediction(size_t sample,
                                                     const Prediction& prediction) const {
  size_t prediction_length = strategy->prediction_length();
  if (prediction.size() != prediction_length) {
    throw std::runtime_error("Prediction for sample " + std::to_string(sample) +
                             " did not have the expected length.");
  }
}

} // namespace grf
