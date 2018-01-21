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
#include "tuning/ParameterTuner.h"
#include "utilities/ForestTestUtilities.h"

#include "catch.hpp"


TEST_CASE("tuning selects a reasonable min node size", "[tuning, forest]") {
  uint outcome_index = 10;
  double alpha = 0.0;

  ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index, alpha);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4, 1);
  ParameterTuner tuner(trainer, predictor, outcome_index);

  Data* data = load_data("test/forest/resources/gaussian_data.csv");
  ForestOptions options = ForestTestUtilities::default_options();

  uint upper_bound = (uint) (data->get_num_rows() * options.get_sample_fraction()) / 4;
  uint min_node_size = tuner.tune_min_node_size(data, options);

  REQUIRE(min_node_size > upper_bound / 2);
  delete data;
}
