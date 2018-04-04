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
#include "forest/ForestTrainers.h"
#include "utilities/ForestTestUtilities.h"

#include "catch.hpp"

TEST_CASE("using a sparse data representation produces the same predictions", "[data]") {
  Data* data = load_data("test/forest/resources/gaussian_data.csv");
  Data* sparse_data = load_sparse_data("test/forest/resources/gaussian_data.csv");

  uint outcome_index = 10;

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 1);
  ForestOptions options = ForestTestUtilities::default_options();

  // Train and predict using the default data format.
  Forest forest = trainer.train(data, options);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  // Train and predict using the sparse data format.
  Forest sparse_forest = trainer.train(sparse_data, options);
  std::vector<Prediction> sparse_predictions = predictor.predict_oob(sparse_forest, sparse_data);

  // Check that the predictions are the same.
  REQUIRE(predictions.size() == sparse_predictions.size());
  for (size_t i = 0; i < predictions.size(); i++) {
    Prediction prediction = predictions[i];
    Prediction sparse_prediction = sparse_predictions[i];

    double value = prediction.get_predictions()[0];
    double sparse_value = sparse_prediction.get_predictions()[0];
    REQUIRE(equal_doubles(value, sparse_value, 1e-5));
  }

  delete data;
  delete sparse_data;
}
