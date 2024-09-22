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

#include "commons/utility.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "utilities/ForestTestUtilities.h"

#include "catch.hpp"

using namespace grf;

TEST_CASE("LLF gives reasonable prediction on friedman data", "[local linear], [forest]") {
  auto data_vec = load_data("test/forest/resources/friedman.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  std::vector<size_t> linear_correction_variables = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> lambda = {0.1};

  bool honesty = true;
  double honesty_fraction = 0.5;
  bool prune = true;
  uint num_trees = 50;
  double sample_fraction = 0.35;
  uint mtry = 3;
  uint min_node_size = 3;
  double alpha = 0.0;
  double imbalance_penalty = 0.0;
  std::vector<size_t> empty_clusters;
  uint samples_per_cluster = 0;
  uint num_threads = 1;
  size_t ci_group_size = 1;
  uint seed = 42;
  ForestOptions options (
      num_trees, ci_group_size, sample_fraction,
      mtry, min_node_size, honesty, honesty_fraction, prune,
      alpha, imbalance_penalty, num_threads, seed, true, empty_clusters, samples_per_cluster);
  ForestTrainer trainer = regression_trainer();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = ll_regression_predictor(
      num_threads, lambda, false, linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, false);

  const std::vector<double>& p = predictions[0].get_predictions();

  REQUIRE(equal_doubles(14.9163, p[0], 3.0));
}

TEST_CASE("LLF predictions vary linearly with Y", "[local linear], [forest]") {
  uint outcome_index = 10;
  auto data_vec = load_data("test/forest/resources/small_gaussian_data.csv");
  Data data(data_vec);
  data.set_outcome_index(outcome_index);

  std::vector<size_t> linear_correction_variables = {1, 4, 7};
  std::vector<double> lambda = {0.1};

  // Run the original forest.
  ForestTrainer trainer = regression_trainer();
  ForestOptions options = ForestTestUtilities::default_honest_options();
  Forest forest = trainer.train(data, options);

  uint num_threads = 1;
  size_t ci_group_size = 1;

  ForestPredictor predictor = ll_regression_predictor(num_threads,
      lambda, false, linear_correction_variables);

  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, false);

  // Shift each outcome by 1, and re-run the forest.
  for (size_t r = 0; r < data.get_num_rows(); r++) {
    double outcome = data.get(r, outcome_index);
    set_data(data_vec, r, outcome_index, outcome + 1);
  }

  Forest shifted_forest = trainer.train(data, options);
  ForestPredictor shifted_predictor = ll_regression_predictor(num_threads,
      lambda, false, linear_correction_variables);
  std::vector<Prediction> shifted_predictions = shifted_predictor.predict_oob(shifted_forest, data, false);

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
}

TEST_CASE("local linear forests give reasonable variance estimates", "[regression, forest]") {
  auto data_vec = load_data("test/forest/resources/gaussian_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);

  double alpha = 0.10;
  double imbalance_penalty = 0.07;

  std::vector<size_t> linear_correction_variables = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> lambda = {0.1};

  bool honesty = true;
  double honesty_fraction = 0.5;
  bool prune = true;
  uint num_trees = 50;
  double sample_fraction = 0.35;
  uint mtry = 3;
  uint min_node_size = 3;
  std::vector<size_t> empty_clusters;
  uint samples_per_cluster = 0;
  uint num_threads = 1;
  size_t ci_group_size = 2;
  uint seed = 42;
  ForestOptions options (
      num_trees, ci_group_size, sample_fraction,
      mtry, min_node_size, honesty, honesty_fraction, prune,
      alpha, imbalance_penalty, num_threads, seed, true, empty_clusters, samples_per_cluster);
  ForestTrainer trainer = regression_trainer();
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = ll_regression_predictor(4, lambda, false, linear_correction_variables);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, true);

  for (const Prediction& prediction : predictions) {
    REQUIRE(prediction.contains_variance_estimates());

    double variance_estimate = prediction.get_variance_estimates()[0];
    REQUIRE(variance_estimate > 0);
  }
}

TEST_CASE("LLF causal predictions are unaffected by shifts in Y", "[local linear], [forest]") {
  auto data_vec = load_data("test/forest/resources/causal_data_ll.csv");
  Data data(data_vec);

  uint outcome_index = 10;
  uint treatment_index = 11;

  data.set_outcome_index(outcome_index);
  data.set_treatment_index(treatment_index);
  data.set_instrument_index(treatment_index);

  std::vector<size_t> linear_correction_variables = {3};
  std::vector<double> lambda = {0.1};

  double reduced_form_weight = 0.0;
  bool stabilize_splits = false;

  ForestTrainer trainer = instrumental_trainer(reduced_form_weight, stabilize_splits);
  ForestOptions options = ForestTestUtilities::default_options();

  Forest forest = trainer.train(data, options);

  uint num_threads = 1;

  ForestPredictor predictor = ll_causal_predictor(num_threads,
      lambda, false, linear_correction_variables);

  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, false);

  // Shift each outcome by 1, and re-run the forest.
  for (size_t r = 0; r < data.get_num_rows(); r++) {
    double outcome = data.get(r, outcome_index);
    set_data(data_vec, r, outcome_index, outcome + 1);
  }

  Forest shifted_forest = trainer.train(data, options);
  ForestPredictor shifted_predictor = ll_causal_predictor(num_threads,
      lambda, false, linear_correction_variables);
  std::vector<Prediction> shifted_predictions = shifted_predictor.predict_oob(shifted_forest, data, false);

  REQUIRE(predictions.size() == shifted_predictions.size());
  double delta = 0.0;
  for (size_t i = 0; i < predictions.size(); i++) {
    Prediction prediction = predictions[i];
    Prediction shifted_prediction = shifted_predictions[i];

    double value = prediction.get_predictions()[0];
    double shifted_value = shifted_prediction.get_predictions()[0];

    delta += shifted_value - value;
  }

  REQUIRE(delta / predictions.size() < 1e-1);
}
