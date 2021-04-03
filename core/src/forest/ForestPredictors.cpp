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
#include "prediction/InstrumentalPredictionStrategy.h"
#include "prediction/MultiCausalPredictionStrategy.h"
#include "prediction/QuantilePredictionStrategy.h"
#include "prediction/ProbabilityPredictionStrategy.h"
#include "prediction/RegressionPredictionStrategy.h"
#include "prediction/MultiRegressionPredictionStrategy.h"
#include "prediction/LocalLinearPredictionStrategy.h"
#include "prediction/LLCausalPredictionStrategy.h"
#include "prediction/SurvivalPredictionStrategy.h"
#include "prediction/CausalSurvivalPredictionStrategy.h"

namespace grf {

ForestPredictor instrumental_predictor(uint num_threads) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor multi_causal_predictor(uint num_threads, size_t num_treatments, size_t num_outcomes) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new MultiCausalPredictionStrategy(num_treatments, num_outcomes));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor quantile_predictor(uint num_threads,
                                   const std::vector<double>& quantiles) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<DefaultPredictionStrategy> prediction_strategy(new QuantilePredictionStrategy(quantiles));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor probability_predictor(uint num_threads, size_t num_classes) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new ProbabilityPredictionStrategy(num_classes));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor regression_predictor(uint num_threads) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor multi_regression_predictor(uint num_threads, size_t num_outcomes) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new MultiRegressionPredictionStrategy(num_outcomes));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor ll_regression_predictor(uint num_threads,
                                        std::vector<double> lambdas,
                                        bool weight_penalty,
                                        std::vector<size_t> linear_correction_variables) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<DefaultPredictionStrategy> prediction_strategy(
      new LocalLinearPredictionStrategy(lambdas, weight_penalty, linear_correction_variables));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor ll_causal_predictor(uint num_threads,
                                    std::vector<double> lambdas,
                                    bool weight_penalty,
                                    std::vector<size_t> linear_correction_variables) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<DefaultPredictionStrategy> prediction_strategy(
          new LLCausalPredictionStrategy(lambdas, weight_penalty, linear_correction_variables));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor survival_predictor(uint num_threads, size_t num_failures, int prediction_type) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<DefaultPredictionStrategy> prediction_strategy(
    new SurvivalPredictionStrategy(num_failures, prediction_type));
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

ForestPredictor causal_survival_predictor(uint num_threads) {
  num_threads = ForestOptions::validate_num_threads(num_threads);
  std::unique_ptr<OptimizedPredictionStrategy> prediction_strategy(new CausalSurvivalPredictionStrategy());
  return ForestPredictor(num_threads, std::move(prediction_strategy));
}

} // namespace grf
