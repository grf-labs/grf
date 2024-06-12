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
#include "prediction/InstrumentalPredictionStrategy.h"
#include "prediction/MultiCausalPredictionStrategy.h"

#include "catch.hpp"

using namespace grf;

TEST_CASE("multi causal predictions with one treatment is identical to causal forest predictions", "[multi_causal, prediction]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(11);
  data.set_instrument_index(11);
  std::vector<std::vector<size_t>> leaf_samples{
    {0, 1, 2, 3, 4, 5},
    {6, 7, 8, 9, 10, 11},
    {12, 13, 14, 15, 16},
    {},
    {21, 22, 38, 41, 18},
    {21, 22, 38},
    {21, 22},
    {87},
    {87, 87},
    {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130}
  };
  size_t num_nodes = leaf_samples.size();

  InstrumentalPredictionStrategy prediction_strategy;
  MultiCausalPredictionStrategy multi_prediction_strategy(1, 1);
  PredictionValues prediction_values = prediction_strategy.precompute_prediction_values(leaf_samples, data);
  PredictionValues multi_prediction_values = multi_prediction_strategy.precompute_prediction_values(leaf_samples, data);

  REQUIRE(prediction_values.get_num_nodes() == multi_prediction_values.get_num_nodes());
  REQUIRE(prediction_values.get_num_nodes() == num_nodes);

  for (size_t i = 0; i < num_nodes; i++) {
    if (prediction_values.empty(i)) {
      REQUIRE(multi_prediction_values.empty(i));
    } else {
      std::vector<double> prediction = prediction_strategy.predict(prediction_values.get_values(i));
      std::vector<double> prediction_multi = multi_prediction_strategy.predict(multi_prediction_values.get_values(i));
      REQUIRE(prediction.size() == prediction_multi.size());
      REQUIRE(equal_doubles(prediction[0], prediction_multi[0], 1e-5));
    }
  }
}

TEST_CASE("multi causal predictions with one continuous treatment is identical to causal forest predictions", "[multi_causal, prediction]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(0); // Set the treatment variable to the first continuous covariate
  data.set_instrument_index(0);

  std::vector<std::vector<size_t>> leaf_samples{
    {0, 1, 2, 3, 4, 5},
    {6, 7, 8, 9, 10, 11},
    {12, 13, 14, 15, 16},
    {},
    {21, 22, 38, 41, 18},
    {21, 22, 38},
    {21, 22},
    {87},
    {87, 87},
    {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130}
  };
  size_t num_nodes = leaf_samples.size();

  InstrumentalPredictionStrategy prediction_strategy;
  MultiCausalPredictionStrategy multi_prediction_strategy(1, 1);
  PredictionValues prediction_values = prediction_strategy.precompute_prediction_values(leaf_samples, data);
  PredictionValues multi_prediction_values = multi_prediction_strategy.precompute_prediction_values(leaf_samples, data);

  REQUIRE(prediction_values.get_num_nodes() == multi_prediction_values.get_num_nodes());
  REQUIRE(prediction_values.get_num_nodes() == num_nodes);

  for (size_t i = 0; i < num_nodes; i++) {
    if (prediction_values.empty(i)) {
      REQUIRE(multi_prediction_values.empty(i));
    } else {
      std::vector<double> prediction = prediction_strategy.predict(prediction_values.get_values(i));
      std::vector<double> prediction_multi = multi_prediction_strategy.predict(multi_prediction_values.get_values(i));
      REQUIRE(prediction.size() == prediction_multi.size());
      REQUIRE(equal_doubles(prediction[0], prediction_multi[0], 1e-5));
    }
  }
}

TEST_CASE("multi causal predictions with one continuous treatment and sample weights is identical to causal forest predictions", "[multi_causal, prediction]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(0); // Set the treatment variable to the first continuous covariate
  data.set_instrument_index(0);
  data.set_weight_index(1); // Use covariate in data column 1 as dummy sample weights

  for(size_t row = 0; row < data.get_num_rows(); row++) {
    double value = data.get(row, 1);
    double weight = value < 0 ? -value : value;
    set_data(data_vec, row, 1, weight);
  }

  std::vector<std::vector<size_t>> leaf_samples{
    {0, 1, 2, 3, 4, 5},
    {6, 7, 8, 9, 10, 11},
    {12, 13, 14, 15, 16},
    {},
    {21, 22, 38, 41, 18},
    {21, 22, 38},
    {21, 22},
    {87},
    {87, 87},
    {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130}
  };
  size_t num_nodes = leaf_samples.size();

  InstrumentalPredictionStrategy prediction_strategy;
  MultiCausalPredictionStrategy multi_prediction_strategy(1, 1);
  PredictionValues prediction_values = prediction_strategy.precompute_prediction_values(leaf_samples, data);
  PredictionValues multi_prediction_values = multi_prediction_strategy.precompute_prediction_values(leaf_samples, data);

  REQUIRE(prediction_values.get_num_nodes() == multi_prediction_values.get_num_nodes());
  REQUIRE(prediction_values.get_num_nodes() == num_nodes);

  for (size_t i = 0; i < num_nodes; i++) {
    if (prediction_values.empty(i)) {
      REQUIRE(multi_prediction_values.empty(i));
    } else {
      std::vector<double> prediction = prediction_strategy.predict(prediction_values.get_values(i));
      std::vector<double> prediction_multi = multi_prediction_strategy.predict(multi_prediction_values.get_values(i));
      REQUIRE(prediction.size() == prediction_multi.size());
      if (std::isinf(prediction[0]) && std::isinf(prediction_multi[0])) {
        continue;
      }
      REQUIRE(equal_doubles(prediction[0], prediction_multi[0], 1e-5));
    }
  }
}

TEST_CASE("multi causal variance estimates with one continuous treatment is identical to causal forest", "[multi_causal, prediction]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(0); // Set the treatment variable to the first continuous covariate
  data.set_instrument_index(0);

  std::vector<std::vector<size_t>> leaf_samples{
    {0, 1, 2, 3, 4, 5},
    {6, 7, 8, 9, 10, 11},
    {12, 13, 14, 15, 16},
    {},
    {21, 22, 38, 41, 18},
    {21, 22, 38},
    {87},
    {87, 87},
    {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130}
  };
  size_t num_nodes = leaf_samples.size();

  InstrumentalPredictionStrategy prediction_strategy;
  MultiCausalPredictionStrategy multi_prediction_strategy(1, 1);
  PredictionValues prediction_values = prediction_strategy.precompute_prediction_values(leaf_samples, data);
  PredictionValues multi_prediction_values = multi_prediction_strategy.precompute_prediction_values(leaf_samples, data);

  REQUIRE(prediction_values.get_num_nodes() == multi_prediction_values.get_num_nodes());
  REQUIRE(prediction_values.get_num_nodes() == num_nodes);

  for (size_t i = 0; i < num_nodes; i++) {
    if (prediction_values.empty(i)) {
      REQUIRE(multi_prediction_values.empty(i));
    } else {
      std::vector<double> prediction = prediction_strategy.compute_variance(prediction_values.get_values(i), prediction_values, 2);
      std::vector<double> prediction_multi = multi_prediction_strategy.compute_variance(multi_prediction_values.get_values(i), multi_prediction_values, 2);
      REQUIRE(prediction.size() == prediction_multi.size());
      REQUIRE(equal_doubles(prediction[0], prediction_multi[0], 1e-5));
    }
  }
}

TEST_CASE("sample weighted multi causal variance estimates with one continuous treatment is identical to causal forest", "[multi_causal, prediction]") {
  auto data_vec = load_data("test/forest/resources/causal_data.csv");
  Data data(data_vec);
  data.set_outcome_index(10);
  data.set_treatment_index(0); // Set the treatment variable to the first continuous covariate
  data.set_instrument_index(0);
  data.set_weight_index(1); // Use covariate in data column 1 as dummy sample weights

  for(size_t row = 0; row < data.get_num_rows(); row++) {
    double value = data.get(row, 1);
    double weight = value < 0 ? -value : value;
    set_data(data_vec, row, 1, weight);
  }

  std::vector<std::vector<size_t>> leaf_samples{
    {0, 1, 2, 3, 4, 5},
    {6, 7, 8, 9, 10, 11},
    {12, 13, 14, 15, 16},
    {},
    {21, 22, 38, 41, 18},
    {21, 22, 38},
    {87},
    {87, 87},
    {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130}
  };
  size_t num_nodes = leaf_samples.size();

  InstrumentalPredictionStrategy prediction_strategy;
  MultiCausalPredictionStrategy multi_prediction_strategy(1, 1);
  PredictionValues prediction_values = prediction_strategy.precompute_prediction_values(leaf_samples, data);
  PredictionValues multi_prediction_values = multi_prediction_strategy.precompute_prediction_values(leaf_samples, data);

  REQUIRE(prediction_values.get_num_nodes() == multi_prediction_values.get_num_nodes());
  REQUIRE(prediction_values.get_num_nodes() == num_nodes);

  for (size_t i = 0; i < num_nodes; i++) {
    if (prediction_values.empty(i)) {
      REQUIRE(multi_prediction_values.empty(i));
    } else {
      std::vector<double> prediction = prediction_strategy.compute_variance(prediction_values.get_values(i), prediction_values, 2);
      std::vector<double> prediction_multi = multi_prediction_strategy.compute_variance(multi_prediction_values.get_values(i), multi_prediction_values, 2);
      REQUIRE(prediction.size() == prediction_multi.size());
      REQUIRE(equal_doubles(prediction[0], prediction_multi[0], 1e-5));
    }
  }
}
