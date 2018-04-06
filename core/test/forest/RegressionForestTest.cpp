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


TEST_CASE("honest regression forests are shift invariant", "[regression, forest]") {
  // Run the original forest.
  Data* data = load_data("test/forest/resources/gaussian_data.csv");
  uint outcome_index = 10;

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index);
  ForestOptions options = ForestTestUtilities::default_honest_options();

  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 1);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  // Shift each outcome by 1, and re-run the forest.
  bool error;
  for (size_t r = 0; r < data->get_num_rows(); r++) {
    double outcome = data->get(r, outcome_index);
    data->set(outcome_index, r, outcome + 1, error);
  }

  Forest shifted_forest = trainer.train(data, options);
  ForestPredictor shifted_predictor = ForestPredictors::regression_predictor(4, 1);
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

TEST_CASE("regression forests give reasonable variance estimates", "[regression, forest]") {
  Data* data = load_data("test/forest/resources/gaussian_data.csv");
  uint outcome_index = 10;
  double alpha = 0.10;
  double imbalance_penalty = 0.07;

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index);
  ForestOptions options = ForestTestUtilities::default_options(false, 2);

  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 2);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  for (size_t i = 0; i < predictions.size(); i++) {
    Prediction prediction = predictions[i];
    REQUIRE(prediction.contains_variance_estimates());

    double variance_estimate = prediction.get_variance_estimates()[0];
    REQUIRE(variance_estimate > 0);
  }

  delete data;
}

TEST_CASE("regression error estimates are shift invariant", "[regression, forest]") {
  // Run the original forest.
  Data* data = load_data("test/forest/resources/gaussian_data.csv");
  uint outcome_index = 10;

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index);
  ForestOptions options = ForestTestUtilities::default_honest_options();

  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 1);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  // Shift each outcome by 1, and re-run the forest.
  bool data_error;
  for (size_t r = 0; r < data->get_num_rows(); r++) {
    double outcome = data->get(r, outcome_index);
    data->set(outcome_index, r, outcome + 1, data_error);
  }

  Forest shifted_forest = trainer.train(data, options);
  ForestPredictor shifted_predictor = ForestPredictors::regression_predictor(4, 1);
  std::vector<Prediction> shifted_predictions = shifted_predictor.predict_oob(shifted_forest, data);

  REQUIRE(predictions.size() == shifted_predictions.size());
  double delta = 0.0;
  for (size_t i = 0; i < predictions.size(); i++) {
    Prediction prediction = predictions[i];
    Prediction shifted_prediction = shifted_predictions[i];

    double error = prediction.get_error_estimates()[0];
    double shifted_error = shifted_prediction.get_error_estimates()[0];

    delta += shifted_error - error;
  }

  REQUIRE(equal_doubles(delta / predictions.size(), 0, 1e-1));
  delete data;
}

