/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include "catch.hpp"

TEST_CASE("forests don't crash when there are fewer trees than threads", "[forest]") {
  uint outcome_index = 10;
  double alpha = 0.10;
  double lambda = 0.07;

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index, alpha, lambda);
  Data* data = load_data("test/forest/resources/gaussian_data.csv");

  uint mtry = 3;
  uint num_trees = 2;
  uint seed = 42;
  uint num_threads = 4;
  uint min_node_size = 1;
  bool sample_with_replacement = true;
  double sample_fraction = 0.35;
  bool honesty = true;
  uint ci_group_size = 2;

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry,
        min_node_size, honesty, sample_with_replacement, num_threads, seed);

  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 2);
  predictor.predict_oob(forest, data);
  delete data;
}
