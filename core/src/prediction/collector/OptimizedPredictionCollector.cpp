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

#include <future>
#include <stdexcept>

#include "prediction/collector/OptimizedPredictionCollector.h"
#include "commons/utility.h"

namespace grf {

OptimizedPredictionCollector::OptimizedPredictionCollector(std::unique_ptr<OptimizedPredictionStrategy> strategy, uint num_threads):
    strategy(std::move(strategy)), num_threads(num_threads) {}

std::vector<Prediction> OptimizedPredictionCollector::collect_predictions(const Forest& forest,
                                                                          const Data& train_data,
                                                                          const Data& data,
                                                                          const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
                                                                          const std::vector<std::vector<bool>>& valid_trees_by_sample,
                                                                          bool estimate_variance,
                                                                          bool estimate_error) const {
  size_t num_samples = data.get_num_rows();
  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, static_cast<uint>(num_samples - 1), num_threads);

  std::vector<std::future<std::vector<Prediction>>> futures;
  futures.reserve(thread_ranges.size());

  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_samples_batch = thread_ranges[i + 1] - start_index;

    futures.push_back(std::async(std::launch::async,
                                 &OptimizedPredictionCollector::collect_predictions_batch,
                                 this,
                                 std::ref(forest),
                                 std::ref(train_data),
                                 std::ref(data),
                                 std::ref(leaf_nodes_by_tree),
                                 std::ref(valid_trees_by_sample),
                                 estimate_variance,
                                 estimate_error,
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

std::vector<Prediction> OptimizedPredictionCollector::collect_predictions_batch(const Forest& forest,
                                                                                const Data& train_data,
                                                                                const Data& data,
                                                                                const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
                                                                                const std::vector<std::vector<bool>>& valid_trees_by_sample,
                                                                                bool estimate_variance,
                                                                                bool estimate_error,
                                                                                size_t start,
                                                                                size_t num_samples) const {
  size_t num_trees = forest.get_trees().size();
  bool record_leaf_values = estimate_variance || estimate_error;

  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  for (size_t sample = start; sample < num_samples + start; ++sample) {
    std::vector<double> average_value;
    std::vector<std::vector<double>> leaf_values;
    if (record_leaf_values) {
      leaf_values.resize(num_trees);
    }

    // Create a list of weighted neighbors for this sample.
    uint num_leaves = 0;
    for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
      if (!valid_trees_by_sample[sample][tree_index]) {
        continue;
      }

      const std::vector<size_t>& leaf_nodes = leaf_nodes_by_tree.at(tree_index);
      size_t node = leaf_nodes.at(sample);

      const std::unique_ptr<Tree>& tree = forest.get_trees()[tree_index];
      const PredictionValues& prediction_values = tree->get_prediction_values();

      if (!prediction_values.empty(node)) {
        num_leaves++;
        add_prediction_values(node, prediction_values, average_value);
        if (record_leaf_values) {
          leaf_values[tree_index] = prediction_values.get_values(node);
        }
      }
    }

    // If this sample has no neighbors, then return placeholder predictions. Note
    // that this can only occur when honesty is enabled, and is expected to be rare.
    if (num_leaves == 0) {
      std::vector<double> nan(strategy->prediction_length(), NAN);
      std::vector<double> nan_error(1, NAN);
      predictions.emplace_back(nan, estimate_variance ? nan : std::vector<double>(), nan_error, nan_error);
      continue;
    }

    normalize_prediction_values(num_leaves, average_value);
    std::vector<double> point_prediction = strategy->predict(average_value);

    PredictionValues prediction_values(leaf_values, strategy->prediction_value_length());
    std::vector<double> variance = estimate_variance
        ? strategy->compute_variance(average_value, prediction_values, forest.get_ci_group_size())
        : std::vector<double>();

    std::vector<double> mse;
    std::vector<double> mce;

    if (estimate_error) {
      std::vector<std::pair<double, double>> error = strategy->compute_error(
              sample, average_value, prediction_values, data);

      mse.push_back(error[0].first);
      mce.push_back(error[0].second);
    }

    Prediction prediction(point_prediction, variance, mse, mce);

    validate_prediction(sample, prediction);
    predictions.push_back(prediction);
  }
  return predictions;
}

void OptimizedPredictionCollector::add_prediction_values(size_t node,
    const PredictionValues& prediction_values,
    std::vector<double>& combined_average) const {
  if (combined_average.empty()) {
    combined_average.resize(prediction_values.get_num_types());
  }

  for (size_t type = 0; type < prediction_values.get_num_types(); ++type) {
    combined_average[type] += prediction_values.get(node, type);
  }
}

void OptimizedPredictionCollector::normalize_prediction_values(size_t num_leaves,
                                                               std::vector<double>& combined_average) const {
  for (double& value : combined_average) {
    value /= num_leaves;
  }
}

void OptimizedPredictionCollector::validate_prediction(size_t sample,
                                                       const Prediction& prediction) const {
  size_t prediction_length = strategy->prediction_length();
  if (prediction.size() != prediction_length) {
    throw std::runtime_error("Prediction for sample " + std::to_string(sample) +
                             " did not have the expected length.");
  }
}

} // namespace grf
