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

#include "commons/utility.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "utilities/ForestTestUtilities.h"

#include "catch.hpp"

using namespace grf;

TEST_CASE("forests don't crash when there are fewer trees than threads", "[forest]") {
  ForestTrainer trainer = regression_trainer();
  auto data_vec = load_data("test/forest/resources/gaussian_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  uint mtry = 3;
  uint num_trees = 2;
  uint seed = 42;
  uint num_threads = 4;
  uint min_node_size = 1;
  size_t ci_group_size = 2;
  double sample_fraction = 0.35;
  bool honesty = true;
  double honesty_fraction = 0.5;
  bool prune = true;
  double alpha = 0.10;
  double imbalance_penalty = 0.07;
  std::vector<size_t> empty_clusters;
  uint samples_per_cluster = 0;

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty, honesty_fraction,
          prune, alpha, imbalance_penalty, num_threads, seed, empty_clusters, samples_per_cluster);

  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = regression_predictor(4);
  predictor.predict_oob(forest, data, true);
}

TEST_CASE("basic forest merges work", "[regression, forest]") {
  auto data_vec = load_data("test/forest/resources/gaussian_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  ForestTrainer trainer = regression_trainer();
  ForestOptions options = ForestTestUtilities::default_options(false, 2);

  std::vector<Forest> forests;
  forests.push_back(trainer.train(data, options));
  forests.push_back(trainer.train(data, options));
  forests.push_back(trainer.train(data, options));

  size_t num_trees = forests[0].get_trees().size();
  size_t num_variables = forests[0].get_num_variables();
  size_t ci_group_size = forests[0].get_ci_group_size();

  Forest big_forest = Forest::merge(forests);

  REQUIRE(num_trees == 50);
  REQUIRE(big_forest.get_trees().size() == 150);

  REQUIRE(num_variables == big_forest.get_num_variables());
  REQUIRE(ci_group_size == big_forest.get_ci_group_size());

  ForestPredictor predictor = regression_predictor(4);
  std::vector<Prediction> predictions = predictor.predict_oob(big_forest, data, false);

  REQUIRE(predictions.size() == data.get_num_rows());
}

TEST_CASE("forests with different ci_group_size cannot be merged", "[regression, forest]") {
  auto data_vec = load_data("test/forest/resources/gaussian_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  ForestTrainer trainer = regression_trainer();
  std::vector<Forest> forests;

  ForestOptions options = ForestTestUtilities::default_options(false, 1);
  forests.push_back(trainer.train(data, options));

  ForestOptions options_with_ci = ForestTestUtilities::default_options(false, 2);
  forests.push_back(trainer.train(data, options_with_ci));

  try {
    Forest big_forest = Forest::merge(forests);
    FAIL();
  } catch (const std::runtime_error&) {
    // Expected exception.
  }
}
