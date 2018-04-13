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

#include "prediction/DefaultPredictionStrategy.h"
#include "commons/utility.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestTrainer.h"
#include "utilities/FileTestUtilities.h"
#include "utilities/ForestTestUtilities.h"

#include "forest/ForestTrainers.h"
#include "forest/ForestPredictors.h"

#include "catch.hpp"

//#define UPDATE_PREDICTION_FILES

bool equal_predictions(const std::vector<Prediction>& actual_predictions,
                       const std::vector<std::vector<double>>& expected_predictions) {
  if (actual_predictions.size() != expected_predictions.size()) {
    return false;
  }

  for (size_t i = 0; i < actual_predictions.size(); ++i) {
    Prediction prediction = actual_predictions[i];
    std::vector<double> expected_prediction = expected_predictions[i];
    if (prediction.size() != expected_prediction.size()) {
      return false;
    }

    for (size_t j = 0; j < prediction.size(); ++j) {
      double value = prediction.get_predictions()[j];
      if (!equal_doubles(value, expected_prediction[j], 1e-2)) {
        return false;
      }
    }
  }

  return true;
}

void update_predictions_file(const std::string file_name, std::vector<Prediction> predictions) {
  std::vector<std::vector<double>> values;
  for (auto& prediction : predictions) {
    values.push_back(prediction.get_predictions());
  }
  FileTestUtilities::write_csv_file(file_name, values);
}

TEST_CASE("quantile forest predictions have not changed", "[quantile], [characterization]") {
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  Data* data = load_data("test/forest/resources/quantile_data.csv");

  ForestTrainer trainer = ForestTrainers::quantile_trainer(10, quantiles);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(4, quantiles);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<Prediction> predictions = predictor.predict(forest, data);


#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/quantile_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/quantile_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/quantile_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/quantile_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));

  delete data;
}

TEST_CASE("causal forest predictions have not changed", "[causal], [characterization]") {
  Data* data = load_data("test/forest/resources/causal_data.csv");

  double reduced_form_weight = 0.0;
  bool stabilize_splits = false;

  ForestTrainer trainer = ForestTrainers::instrumental_trainer(
      10, 11, 11, reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(4, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<Prediction> predictions = predictor.predict(forest, data);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/causal_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/causal_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));

  delete data;
}

TEST_CASE("causal forest predictions with stable splitting have not changed", "[causal], [characterization]") {
  Data* data = load_data("test/forest/resources/causal_data.csv");

  double reduced_form_weight = 0.0;
  bool stabilize_splits = true;

  ForestTrainer trainer = ForestTrainers::instrumental_trainer(
      10, 11, 11, reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(4, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<Prediction> predictions = predictor.predict(forest, data);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/stable_causal_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/stable_causal_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_causal_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_causal_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));

  delete data;
}

TEST_CASE("regression forest predictions have not changed", "[regression], [characterization]") {
  Data* data = load_data("test/forest/resources/regression_data.csv");

  ForestTrainer trainer = ForestTrainers::regression_trainer(10);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<Prediction> predictions = predictor.predict(forest, data);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/regression_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/regression_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));

  delete data;
}
