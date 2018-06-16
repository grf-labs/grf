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

#include "commons/utility.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "utilities/ForestTestUtilities.h"
#include <chrono>
#include <prediction/LocalLinearPredictionStrategy.h>
#include <prediction/collector/DefaultPredictionCollector.h>

#include "catch.hpp"

TEST_CASE("LLFBenchmark") {
  Data* data = load_data("test/forest/resources/friedman.csv");
  uint outcome_index = 10;
  std::vector<size_t> linear_correction_variables = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index);
  bool honesty = true;
  uint num_trees = 50;
  double sample_fraction = 0.35;
  uint mtry = 3;
  uint min_node_size = 3;
  double alpha = 0.0;
  double imbalance_penalty = 0.0;
  std::vector<size_t> empty_clusters;
  uint samples_per_cluster = 0;
  uint num_threads = 1;
  uint seed = 42;
  ForestOptions options (
      num_trees, 1, sample_fraction, mtry, min_node_size, honesty,
      alpha, imbalance_penalty, num_threads, seed, empty_clusters, samples_per_cluster
  );
  // Benchmark Flags
  size_t num_iter = 1;
  bool track_statistics = false;
  bool print_runtime = false;


  auto start = std::chrono::high_resolution_clock::now();
  Forest forest = trainer.train(data, options);
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;

  Data* queries = data;
  auto strategy = std::make_shared<LocalLinearPredictionStrategy>(
      data, queries,
      .1, false,
      linear_correction_variables
  );
  DefaultPredictionCollector prediction_collector (strategy);
  TreeTraverser tree_traverser (1);
  std::vector<Prediction> predictions;
  size_t num_queries = queries->get_num_rows();
  double avg_leaves, avg_samples_leaf, avg_samples_query;

  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < num_iter; i++) {
    std::vector<std::vector<size_t>> leaf_nodes_by_tree = tree_traverser.get_leaf_nodes(forest, data, queries);
    std::vector<std::vector<bool>> trees_by_sample = tree_traverser.get_valid_trees_by_sample(forest, data, queries);

    if (track_statistics) {
      size_t sum_leaves = 0;
      size_t sum_samples = 0;
      size_t sum_distinct_samples = 0;
      for (size_t sample = 0; sample < num_queries; sample++) {
        std::set<size_t> query_samples;
        for (size_t tree_index = 0; tree_index < forest.get_trees().size(); ++tree_index) {
          if (!trees_by_sample[sample][tree_index]) {
            continue;
          }
          sum_leaves++;

          const std::vector<size_t> &leaf_nodes = leaf_nodes_by_tree.at(tree_index);
          size_t node = leaf_nodes.at(sample);

          std::shared_ptr<Tree> tree = forest.get_trees()[tree_index];
          const std::vector<size_t> &samples = tree->get_leaf_samples()[node];
          sum_samples += samples.size();
          query_samples.insert(samples.begin(), samples.end());
        }
        sum_distinct_samples += query_samples.size();
      }
      avg_leaves = sum_leaves * 1.0 / num_queries;
      avg_samples_leaf = sum_samples * 1.0 / sum_leaves;
      avg_samples_query = sum_distinct_samples * 1.0 / num_queries;
    }

    predictions = prediction_collector.collect_predictions(
        forest, data,
        leaf_nodes_by_tree, trees_by_sample,
        queries
    );
  }
  finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  if (print_runtime) {
    std::cout << "Predicted in: " << (elapsed.count()) / num_iter << std::endl;
    std::cout << "Avg Leafs per Query: " << avg_leaves << std::endl;
    std::cout << "Avg Samples per Leaf: " << avg_samples_leaf << std::endl;
    std::cout << "Avg Samples per Query: " << avg_samples_query << std::endl;
  }

  const std::vector<double>& p = predictions[0].get_predictions();
  REQUIRE(equal_doubles(14.9163, p[0], 1e-3));
}


TEST_CASE("LLF predictions vary linearly with Y", "[local_linear, forest]") {
  Data* data = load_data("test/forest/resources/small_gaussian_data.csv");
  uint outcome_index = 10;
  std::vector<size_t> linear_correction_variables = {1, 4, 7};

  // Run the original forest.
  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index);
  ForestOptions options = ForestTestUtilities::default_honest_options();
  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = ForestPredictors::local_linear_predictor(4, data, data, 0.1, false, linear_correction_variables);
  
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);
  const std::vector<double>& p = predictions[0].get_predictions();
  REQUIRE(equal_doubles(0.0107178, p[0], 1e-5));

  // Shift each outcome by 1, and re-run the forest.
  bool error;
  for (size_t r = 0; r < data->get_num_rows(); r++) {
    double outcome = data->get(r, outcome_index);
    data->set(outcome_index, r, outcome + 1, error);
  }
  
  Forest shifted_forest = trainer.train(data, options);
  ForestPredictor shifted_predictor = ForestPredictors::local_linear_predictor(4, data, data, 0.1, false, linear_correction_variables);
  std::vector<Prediction> shifted_predictions = shifted_predictor.predict_oob(shifted_forest, data);
  
  REQUIRE(predictions.size() == shifted_predictions.size());
  double delta = 0.0;
  for (size_t i = 0; i < predictions.size(); i++) {
    Prediction prediction = predictions[i];
    Prediction shifted_prediction = shifted_predictions[i];
  
    double value = prediction.get_predictions()[0];
    double shifted_value = shifted_prediction.get_predictions()[0];
  
    delta += shifted_value - value;
  }
  
  REQUIRE(equal_doubles(delta / predictions.size(), 1, 1e-1));
  delete data;
}
