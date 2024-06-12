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

TEST_CASE("causal forests are invariant to rescaling of the sample weights", "[causal, forest]") {
  // Run the original forest.
  // we'll overwrite a covariate in the original data with sample weights so we needn't resize the data.
  size_t weight_index = 9;
  size_t outcome_index = 10;
  size_t treatment_index = 11;
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_weight_index(weight_index);
  data.set_outcome_index(outcome_index);
  data.set_treatment_index(treatment_index);
  data.set_instrument_index(treatment_index);

  for(size_t r = 0; r < data.get_num_rows(); r++) {
    double weight = 1.0 / (1.0 + exp(- data.get(r, 1)));
    set_data(data_vec, r, weight_index, weight);
  }

  ForestTrainer trainer = instrumental_trainer(0, true);
  ForestOptions options = ForestTestUtilities::default_honest_options();

  Forest forest = trainer.train(data, options);
  ForestPredictor predictor = instrumental_predictor(4);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, false);

  // Scale weights by n and re-run the forest.
  for (size_t r = 0; r < data.get_num_rows(); r++) {
    double weight = data.get_weight(r) * data.get_num_rows();
    set_data(data_vec, r, weight_index, weight);
  }

  Forest shifted_forest = trainer.train(data, options);
  ForestPredictor shifted_predictor = instrumental_predictor(4);
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

  REQUIRE(equal_doubles(delta / predictions.size(), 0, 1e-1));
}
