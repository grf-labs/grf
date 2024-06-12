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

#include "commons/Data.h"
#include "commons/utility.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "prediction/MultiRegressionPredictionStrategy.h"

#include "catch.hpp"

using namespace grf;

TEST_CASE("multi regression predictions with one outcome is identical to regression predictions", "[multi_regression, prediction]") {
  auto data_vec = load_data("test/forest/resources/regression_data_MIA.csv");
  Data data(data_vec);
  data.set_outcome_index(5);
  std::vector<std::vector<size_t>> leaf_samples{
    {0, 1, 2, 3, 4, 5},
    {6, 7, 8, 9, 10, 11},
    {12, 13, 14, 15, 16},
    {},
    {21, 22, 38, 41, 18},
    {87}
  };
  size_t num_nodes = leaf_samples.size();

  RegressionPredictionStrategy prediction_strategy;
  MultiRegressionPredictionStrategy multi_prediction_strategy(1);
  PredictionValues reg_prediction_values = prediction_strategy.precompute_prediction_values(leaf_samples, data);
  PredictionValues multi_reg_prediction_values = multi_prediction_strategy.precompute_prediction_values(leaf_samples, data);

  REQUIRE(reg_prediction_values.get_num_nodes() == multi_reg_prediction_values.get_num_nodes());
  REQUIRE(reg_prediction_values.get_num_types() == multi_reg_prediction_values.get_num_types());
  REQUIRE(reg_prediction_values.get_num_nodes() == num_nodes);
  for (size_t i = 0; i < num_nodes; i++) {
    if (reg_prediction_values.empty(i)) {
      REQUIRE(multi_reg_prediction_values.empty(i));
    } else {
      std::vector<double> prediction = prediction_strategy.predict(reg_prediction_values.get_values(i));
      std::vector<double> prediction_multi = multi_prediction_strategy.predict(multi_reg_prediction_values.get_values(i));
      REQUIRE(prediction.size() == prediction_multi.size());
      REQUIRE(equal_doubles(prediction[0], prediction_multi[0], 1e-10));
    }
  }
}
