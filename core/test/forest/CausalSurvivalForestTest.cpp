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

TEST_CASE("causal survival forests give positive variance estimates", "[causal survival]") {
  auto data_vec = load_data("test/forest/resources/causal_survival_data.csv");
  Data data(data_vec);
  data.set_treatment_index(5);
  data.set_instrument_index(5);
  data.set_censor_index(6);
  data.set_causal_survival_numerator_index(7);
  data.set_causal_survival_denominator_index(8);

  ForestTrainer trainer = causal_survival_trainer(true);
  ForestOptions options = ForestTestUtilities::default_options(true, 2);
  Forest forest = trainer.train(data, options);

  ForestPredictor predictor = causal_survival_predictor(4);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, true);

  for (const Prediction& prediction : predictions) {
    REQUIRE(prediction.contains_variance_estimates());

    double variance_estimate = prediction.get_variance_estimates()[0];
    REQUIRE(variance_estimate > 0);
  }
}
