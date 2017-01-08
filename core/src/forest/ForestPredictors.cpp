#include "ForestPredictors.h"
#include "InstrumentalPredictionStrategy.h"
#include "QuantilePredictionStrategy.h"
#include "RegressionPredictionStrategy.h"

ForestPredictor ForestPredictors::instrumental_predictor(uint num_threads) {
  std::shared_ptr<PredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());
  return ForestPredictor(num_threads, prediction_strategy);
}

ForestPredictor ForestPredictors::quantile_predictor(uint num_threads,
                                                     const std::vector<double> &quantiles) {
  std::shared_ptr<PredictionStrategy> prediction_strategy(new QuantilePredictionStrategy(quantiles));
  return ForestPredictor(num_threads, prediction_strategy);
}

ForestPredictor ForestPredictors::regression_predictor(uint num_threads) {
  std::shared_ptr<PredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());
  return ForestPredictor(num_threads, prediction_strategy);
}
