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

#include "commons/Data.h"
#include "commons/utility.h"
#include "prediction/SurvivalPredictionStrategy.h"

#include "catch.hpp"

using namespace grf;

TEST_CASE("Kaplan-Meier survival estimates are correct", "[survival], [prediction]") {
  size_t num_failures = 24;
  size_t num_rows = 50;
  size_t num_cols = 2;
  size_t outcome_index = 0;

  std::vector<double> data_matrix = {
    10L, 22L, 19L, 0L, 18L, 7L, 6L, 13L, 4L, 14L, 5L, 10L, 24L,
    4L, 9L, 23L, 4L, 3L, 16L, 11L, 11L, 7L, 20L, 7L, 21L, 1L, 23L,
    10L, 24L, 7L, 15L, 2L, 12L, 8L, 17L, 14L, 9L, 10L, 2L, 11L, 23L,
    20L, 16L, 8L, 8L, 10L, 24L, 23L, 22L, 10L, 0L, 1L, 1L, 0L, 1L,
    1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L,
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L,
    1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L
  };

  std::vector<double> expected_predictions = {
    0.979591836734694, 0.959183673469388, 0.938331854480923, 0.917480035492458,
    0.895635272742637, 0.873790509992817, 0.851945747242997, 0.828280587597358,
    0.803181175851983, 0.77727210566321, 0.746181221436681, 0.712263893189559,
    0.678346564942438, 0.644429236695316, 0.608627612434465, 0.572825988173614,
    0.53463758896204, 0.496449189750465, 0.458260790538891, 0.420072391327317,
    0.378065152194585, 0.336057913061853, 0.288049639767303, 0.192033093178202
  };

  DefaultData data(data_matrix, num_rows, num_cols);
  data.set_outcome_index(outcome_index);
  data.set_censor_index(outcome_index + 1);

  std::unordered_map<size_t, double> weights_by_sample;
  for (size_t i = 0; i < num_rows; i++) {
    weights_by_sample[i] = 1.0;
  }

  SurvivalPredictionStrategy prediction_strategy(num_failures);
  std::vector<double> predictions = prediction_strategy.predict(0, weights_by_sample, data, data);

  for (size_t i = 0; i < predictions.size(); i++) {
    REQUIRE(equal_doubles(predictions[i], expected_predictions[i], 1e-10));
  }
}

TEST_CASE("Kaplan-Meier estimates on duplicated data is the same as with sample weights equal to two", "[survival], [prediction]") {
  size_t num_failures = 24;
  size_t num_rows = 50;
  size_t num_cols = 2;
  size_t outcome_index = 0;

  std::vector<double> data_matrix = {
    10L, 22L, 19L, 0L, 18L, 7L, 6L, 13L, 4L, 14L, 5L, 10L, 24L,
    4L, 9L, 23L, 4L, 3L, 16L, 11L, 11L, 7L, 20L, 7L, 21L, 1L, 23L,
    10L, 24L, 7L, 15L, 2L, 12L, 8L, 17L, 14L, 9L, 10L, 2L, 11L, 23L,
    20L, 16L, 8L, 8L, 10L, 24L, 23L, 22L, 10L, 0L, 1L, 1L, 0L, 1L,
    1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L,
    0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L,
    1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L
  };

  // Duplicate the first 10 samples
  std::vector<double> data_matrix_duplicated = {
    10L, 22L, 19L, 0L, 18L, 7L, 6L, 13L, 4L, 14L, 5L, 10L, 24L,
    4L, 9L, 23L, 4L, 3L, 16L, 11L, 11L, 7L, 20L, 7L, 21L, 1L, 23L,
    10L, 24L, 7L, 15L, 2L, 12L, 8L, 17L, 14L, 9L, 10L, 2L, 11L, 23L,
    20L, 16L, 8L, 8L, 10L, 24L, 23L, 22L, 10L,
    10L, 22L, 19L, 0L, 18L, 7L, 6L, 13L, 4L, 14L, // duplicated
    0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L,
    0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L,
    1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L,
    0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 0L // duplicated
  };

  // add sample weights
  size_t num_duplicates = 10;
  for (size_t i = 0; i < num_rows; i++) {
    if (i < num_duplicates) {
      data_matrix.push_back(2);
    } else {
      data_matrix.push_back(1);
    }
  }

  DefaultData data(data_matrix, num_rows, num_cols + 1);
  data.set_outcome_index(outcome_index);
  data.set_censor_index(outcome_index + 1);
  data.set_weight_index(outcome_index + 2);

  DefaultData data_duplicated(data_matrix_duplicated, num_rows + num_duplicates, num_cols);
  data_duplicated.set_outcome_index(outcome_index);
  data_duplicated.set_censor_index(outcome_index + 1);

  std::unordered_map<size_t, double> weights_by_sample;
  for (size_t i = 0; i < num_rows; i++) {
    weights_by_sample[i] = 1.0;
  }

  SurvivalPredictionStrategy prediction_strategy(num_failures);
  std::vector<double> predictions_weighted = prediction_strategy.predict(0, weights_by_sample, data, data);
  for (size_t i = num_rows; i < num_rows + num_duplicates; i++) {
    weights_by_sample[i] = 1.0;
  }
  std::vector<double> predictions_duplicated = prediction_strategy.predict(0, weights_by_sample, data_duplicated, data_duplicated);

  for (size_t i = 0; i < predictions_duplicated.size(); i++) {
    REQUIRE(equal_doubles(predictions_weighted[i], predictions_duplicated[i], 1e-10));
  }
}
