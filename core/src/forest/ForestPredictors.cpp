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

#include "forest/ForestPredictors.h"
#include "prediction/CustomPredictionStrategy.h"
#include "prediction/InstrumentalPredictionStrategy.h"
#include "prediction/QuantilePredictionStrategy.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "prediction/LocalLinearPredictionStrategy.h"
#include "ForestOptions.h"

ForestPredictor ForestPredictors::custom_predictor(uint num_threads) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::shared_ptr<DefaultPredictionStrategy> prediction_strategy(new CustomPredictionStrategy());
  return ForestPredictor(num_threads, prediction_strategy);
}

ForestPredictor ForestPredictors::instrumental_predictor(uint num_threads,
                                                         uint ci_group_size) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());
  return ForestPredictor(num_threads, ci_group_size, prediction_strategy);
}

ForestPredictor ForestPredictors::quantile_predictor(uint num_threads,
                                                     const std::vector<double>& quantiles) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::shared_ptr<DefaultPredictionStrategy> prediction_strategy(new QuantilePredictionStrategy(quantiles));
  return ForestPredictor(num_threads, prediction_strategy);
}

ForestPredictor ForestPredictors::regression_predictor(uint num_threads,
                                                       uint ci_group_size) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::shared_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());
  return ForestPredictor(num_threads, ci_group_size, prediction_strategy);
}

ForestPredictor ForestPredictors::local_linear_predictor(uint num_threads,
                                                         const Data*original_data,
                                                         const Data *test_data,
                                                         std::vector<double> lambdas,
                                                         bool use_unweighted_penalty,
                                                         std::vector<size_t> linear_correction_variables) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::shared_ptr<DefaultPredictionStrategy> prediction_strategy(new LocalLinearPredictionStrategy(original_data,
                                                                                                   test_data,
                                                                                                   lambdas,
                                                                                                   use_unweighted_penalty,
                                                                                                   linear_correction_variables));
  return ForestPredictor(num_threads, prediction_strategy);
}
