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

using namespace grf;

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

void update_predictions_file(const std::string& file_name,
                             const std::vector<Prediction>& predictions) {
  std::vector<std::vector<double>> values;
  values.reserve(predictions.size());
  for (const auto& prediction : predictions) {
    values.push_back(prediction.get_predictions());
  }
  FileTestUtilities::write_csv_file(file_name, values);
}

TEST_CASE("quantile forest predictions have not changed", "[quantile], [characterization]") {
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  auto data_vec = load_data("test/forest/resources/quantile_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  ForestTrainer trainer = quantile_trainer(quantiles);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = quantile_predictor(4, quantiles);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

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
}

TEST_CASE("quantile forest predictions with NaNs have not changed", "[NaN], [quantile], [characterization]") {
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  auto data_vec = load_data("test/forest/resources/quantile_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  ForestTrainer trainer = quantile_trainer(quantiles);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = quantile_predictor(4, quantiles);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/quantile_oob_predictions_MIA.csv", oob_predictions);
  update_predictions_file("test/forest/resources/quantile_predictions_MIA.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/quantile_oob_predictions_MIA.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/quantile_predictions_MIA.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("probability forest predictions have not changed", "[probability], [characterization]") {
  auto data_vec = load_data("test/forest/resources/probability_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  size_t num_classes = 6;

  ForestTrainer trainer = probability_trainer(num_classes);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = probability_predictor(4, num_classes);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/probability_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/probability_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/probability_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/probability_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}


TEST_CASE("causal forest predictions have not changed", "[causal], [characterization]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);

  double reduced_form_weight = 0.0;
  bool stabilize_splits = false;

  ForestTrainer trainer = instrumental_trainer(reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = instrumental_predictor(4);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

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
}

TEST_CASE("causal forest predictions with stable splitting have not changed", "[causal], [characterization]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);

  double reduced_form_weight = 0.0;
  bool stabilize_splits = true;

  ForestTrainer trainer = instrumental_trainer(reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = instrumental_predictor(4);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

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
}

TEST_CASE("causal forest predictions with sample weights and stable splitting have not changed", "[causal], [characterization]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  size_t weight_index = 9;
  data.set_weight_index(weight_index);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);

  // Use covariate in data column 9 as dummy sample weights
  for(size_t r = 0; r < data.get_num_rows(); r++) {
    double value = data.get(r, weight_index);
    double weight = value < 0 ? -value : value;
    set_data(data_vec, r, weight_index, weight);
  }

  double reduced_form_weight = 0.0;
  bool stabilize_splits = true;

  ForestTrainer trainer = instrumental_trainer(reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = instrumental_predictor(4);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/causal_oob_predictions_sample_weights.csv", oob_predictions);
  update_predictions_file("test/forest/resources/causal_predictions_sample_weights.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_oob_predictions_sample_weights.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_predictions_sample_weights.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("causal forest predictions with NaNs and stable splitting have not changed", "[NaN], [causal], [characterization]") {
  auto data_vec = load_data("test/forest/resources/causal_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);

  double reduced_form_weight = 0.0;
  bool stabilize_splits = true;

  ForestTrainer trainer = instrumental_trainer(reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = instrumental_predictor(4);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/stable_causal_oob_predictions_MIA.csv", oob_predictions);
  update_predictions_file("test/forest/resources/stable_causal_predictions_MIA.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_causal_oob_predictions_MIA.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_causal_predictions_MIA.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("multi causal forest predictions with sample weights have not changed", "[multi causal], [characterization]") {
  auto data_vec = load_data("test/forest/resources/multi_causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(5);
  data.set_treatment_index({6, 7});
  data.set_weight_index(8);

  size_t num_treatments = 2;
  ForestTrainer trainer = multi_causal_trainer(num_treatments, 1, false);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = multi_causal_predictor(4, num_treatments, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/stable_multi_causal_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/stable_multi_causal_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_multi_causal_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_multi_causal_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("multi causal forest predictions with sample weights and stable splitting have not changed", "[multi causal], [characterization]") {
  auto data_vec = load_data("test/forest/resources/multi_causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(5);
  data.set_treatment_index({6, 7});
  data.set_weight_index(8);

  size_t num_treatments = 2;
  ForestTrainer trainer = multi_causal_trainer(num_treatments, 1, true);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = multi_causal_predictor(4, num_treatments, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/stable_multi_causal_oob_predictions_stable_split.csv", oob_predictions);
  update_predictions_file("test/forest/resources/stable_multi_causal_predictions_stable_split.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_multi_causal_oob_predictions_stable_split.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/stable_multi_causal_predictions_stable_split.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("regression forest predictions have not changed", "[regression], [characterization]") {
  auto data_vec = load_data("test/forest/resources/regression_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  ForestTrainer trainer = regression_trainer();
  ForestTrainer multi_trainer = multi_regression_trainer(1);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);
  Forest multi_forest = multi_trainer.train(data, options);

  ForestPredictor predictor = regression_predictor(4);
  ForestPredictor multi_predictor = multi_regression_predictor(4, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);
  std::vector<Prediction> multi_oob_predictions = multi_predictor.predict_oob(multi_forest, data, false);
  std::vector<Prediction> multi_predictions = multi_predictor.predict(multi_forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/regression_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/regression_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));
  REQUIRE(equal_predictions(multi_oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
  REQUIRE(equal_predictions(multi_predictions, expected_predictions));
}

TEST_CASE("regression forest predictions with sample weights have not changed", "[regression], [characterization]") {
  auto data_vec = load_data("test/forest/resources/regression_data.csv");
  Data data(data_vec);
  size_t weight_index = 9;
  data.set_weight_index(weight_index);
  data.set_outcome_index(10);

  // Use covariate in data column 9 as dummy sample weights
  for(size_t r = 0; r < data.get_num_rows(); r++) {
    double value = data.get(r, weight_index);
    double weight = value < 0 ? -value : value;
    set_data(data_vec, r, weight_index, weight);
  }

  ForestTrainer trainer = regression_trainer();
  ForestTrainer multi_trainer = multi_regression_trainer(1);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);
  Forest multi_forest = multi_trainer.train(data, options);

  ForestPredictor predictor = regression_predictor(4);
  ForestPredictor multi_predictor = multi_regression_predictor(4, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);
  std::vector<Prediction> multi_oob_predictions = multi_predictor.predict_oob(multi_forest, data, false);
  std::vector<Prediction> multi_predictions = multi_predictor.predict(multi_forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/regression_oob_predictions_sample_weights.csv", oob_predictions);
  update_predictions_file("test/forest/resources/regression_predictions_sample_weights.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_oob_predictions_sample_weights.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));
  REQUIRE(equal_predictions(multi_oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_predictions_sample_weights.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
  REQUIRE(equal_predictions(multi_predictions, expected_predictions));
}

TEST_CASE("regression forest predictions with NaNs have not changed", "[NaN], [regression], [characterization]") {
  auto data_vec = load_data("test/forest/resources/regression_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(5);

  ForestTrainer trainer = regression_trainer();
  ForestTrainer multi_trainer = multi_regression_trainer(1);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);
  Forest multi_forest = multi_trainer.train(data, options);

  ForestPredictor predictor = regression_predictor(4);
  ForestPredictor multi_predictor = multi_regression_predictor(4, 1);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);
  std::vector<Prediction> multi_oob_predictions = multi_predictor.predict_oob(multi_forest, data, false);
  std::vector<Prediction> multi_predictions = multi_predictor.predict(multi_forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/regression_oob_predictions_MIA.csv", oob_predictions);
  update_predictions_file("test/forest/resources/regression_predictions_MIA.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_oob_predictions_MIA.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));
  REQUIRE(equal_predictions(multi_oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/regression_predictions_MIA.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
  REQUIRE(equal_predictions(multi_predictions, expected_predictions));
}

TEST_CASE("local linear regression forest predictions have not changed",
          "[local linear], [regression], [characterization]") {
  auto data_vec = load_data("test/forest/resources/regression_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  ForestTrainer trainer = regression_trainer();
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  std::vector<double> lambdas = {0, 0.1, 1, 10, 100};
  bool weight_penalty = false;
  std::vector<size_t> linear_correction_variables = {1, 3, 5};
  ForestPredictor predictor = ll_regression_predictor(
      4, lambdas, weight_penalty, linear_correction_variables);

  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data,  false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/ll_regression_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/ll_regression_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/ll_regression_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/ll_regression_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("survival forest predictions have not changed", "[survival], [characterization]") {
  size_t num_failures = 149;
  auto data_vec = load_data("test/forest/resources/survival_data.csv");
  Data data(data_vec);
  data.set_outcome_index(5);
  data.set_censor_index(6);

  ForestTrainer trainer = survival_trainer();
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  int prediction_type = 0;
  ForestPredictor predictor = survival_predictor(4, num_failures, prediction_type);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/survival_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/survival_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/survival_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/survival_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("survival forest predictions with NaNs have not changed", "[NaN], [survival], [characterization]") {
  size_t num_failures = 149;
  auto data_vec = load_data("test/forest/resources/survival_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(5);
  data.set_censor_index(6);

  ForestTrainer trainer = survival_trainer();
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  int prediction_type = 0;
  ForestPredictor predictor = survival_predictor(4, num_failures, prediction_type);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/survival_oob_predictions_MIA.csv", oob_predictions);
  update_predictions_file("test/forest/resources/survival_predictions_MIA.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/survival_oob_predictions_MIA.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/survival_predictions_MIA.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("causal survival forest predictions have not changed", "[causal survival], [characterization]") {
  auto data_vec = load_data("test/forest/resources/causal_survival_data.csv");
  Data data(data_vec);
  data.set_treatment_index(5);
  data.set_instrument_index(5);
  data.set_censor_index(6);
  data.set_causal_survival_numerator_index(7);
  data.set_causal_survival_denominator_index(8);

  ForestTrainer trainer = causal_survival_trainer(true);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = causal_survival_predictor(4);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/causal_survival_oob_predictions.csv", oob_predictions);
  update_predictions_file("test/forest/resources/causal_survival_predictions.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_survival_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_survival_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("causal survival forest predictions with NaNs have not changed", "[NaN], [causal survival], [characterization]") {
  auto data_vec = load_data("test/forest/resources/causal_survival_data_MIA.csv");
  Data data(data_vec);
  data.set_treatment_index(5);
  data.set_instrument_index(5);
  data.set_censor_index(6);
  data.set_causal_survival_numerator_index(7);
  data.set_causal_survival_denominator_index(8);

  ForestTrainer trainer = causal_survival_trainer(true);
  ForestOptions options = ForestTestUtilities::default_options();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = causal_survival_predictor(4);
  std::vector<Prediction> oob_predictions = predictor.predict_oob(forest, data, false);
  std::vector<Prediction> predictions = predictor.predict(forest, data, data, false);

#ifdef UPDATE_PREDICTION_FILES
  update_predictions_file("test/forest/resources/causal_survival_oob_predictions_MIA.csv", oob_predictions);
  update_predictions_file("test/forest/resources/causal_survival_predictions_MIA.csv", predictions);
#endif

  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_survival_oob_predictions_MIA.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::read_csv_file(
      "test/forest/resources/causal_survival_predictions_MIA.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}
