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

#include "prediction/collector/DefaultPredictionCollector.h"

namespace grf {

DefaultPredictionCollector::DefaultPredictionCollector(std::unique_ptr<DefaultPredictionStrategy> strategy):
    strategy(std::move(strategy)) {}

std::vector<Prediction> DefaultPredictionCollector::collect_predictions(
    const Forest& forest,
    const Data& train_data,
    const Data& data,
    const std::vector<std::vector<size_t>>& leaf_nodes_by_tree,
    const std::vector<std::vector<bool>>& valid_trees_by_sample,
    bool estimate_variance,
    bool estimate_error) const {

  size_t num_samples = data.get_num_rows();
  std::vector<Prediction> predictions;
  predictions.reserve(num_samples);

  size_t num_trees = forest.get_trees().size();
  bool record_leaf_samples = estimate_variance || estimate_error;

  for (size_t sample = 0; sample < num_samples; sample++) {
    std::unordered_map<size_t, double> weights_by_sample = weight_computer.compute_weights(
        sample, forest, leaf_nodes_by_tree, valid_trees_by_sample);
    std::vector<std::vector<size_t>> samples_by_tree;

    if (record_leaf_samples) {
      samples_by_tree.resize(num_trees);

      for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
        if (!valid_trees_by_sample[sample][tree_index]) {
          continue;
        }
        const std::vector<size_t>& leaf_nodes = leaf_nodes_by_tree.at(tree_index);
        size_t node = leaf_nodes.at(sample);

        const std::unique_ptr<Tree>& tree = forest.get_trees()[tree_index];
        std::vector<std::vector<size_t>> leaf_samples = tree->get_leaf_samples();
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
