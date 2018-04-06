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

#include "catch.hpp"


TEST_CASE("custom forests predict 0 by default", "[custom, forest]") {
  // Train an honest custom forest.
  Data* data = load_data("test/forest/resources/gaussian_data.csv");
  uint outcome_index = 10;

  ForestTrainer trainer = ForestTrainers::custom_trainer(outcome_index);
  ForestOptions options = ForestTestUtilities::default_honest_options();
  Forest forest = trainer.train(data, options);

  // Predict on the same data.
  ForestPredictor predictor = ForestPredictors::custom_predictor(4);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  // Check the dummy predictions look as expected.
  REQUIRE(predictions.size() == data->get_num_rows());
  for (Prediction prediction : predictions) {
    double value = prediction.get_predictions()[0];
    REQUIRE(equal_doubles(value, 0.0, 1e-10));
  }

  delete data;
}
