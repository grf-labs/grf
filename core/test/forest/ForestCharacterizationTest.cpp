#include "PredictionStrategy.h"
#include "utility.h"
#include "ForestPredictor.h"
#include "ForestTrainer.h"
#include "QuantileRelabelingStrategy.h"
#include "QuantilePredictionStrategy.h"
#include "RegressionSplittingRuleFactory.h"
#include "ProbabilitySplittingRuleFactory.h"
#include "InstrumentalRelabelingStrategy.h"
#include "NoopRelabelingStrategy.h"
#include "InstrumentalPredictionStrategy.h"
#include "RegressionPredictionStrategy.h"
#include "FileTestUtilities.h"
#include "ForestTrainers.h"

#include "catch.hpp"

void init_trainer(ForestTrainer& trainer) {
  uint mtry = 3;
  uint num_trees = 4;
  std::ostream* verbose_out = &std::cout;
  uint seed = 42;
  uint num_threads = 4;
  std::string load_forest_filename = "";
  uint min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  std::vector<size_t> no_split_variables;
  std::string split_select_weights_file = "";
  std::vector<std::string> always_split_variable_names;
  bool sample_with_replacement = true;
  bool memory_saving_splitting = false;
  std::string case_weights_file = "";
  double sample_fraction = 1;

  trainer.init(mtry, num_trees, verbose_out, seed, num_threads, load_forest_filename,
                min_node_size, no_split_variables, split_select_weights_file, always_split_variable_names,
                sample_with_replacement, memory_saving_splitting, case_weights_file, sample_fraction);
}

ForestPredictor create_forest_predictor(std::shared_ptr<PredictionStrategy> prediction_strategy) {
  uint num_threads = 4;
  return ForestPredictor(num_threads, prediction_strategy);
}

TEST_CASE("quantile forest predictions have not changed", "[quantile, characterization]") {
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  Data *data = loadDataFromFile("test/forest/resources/quantile_test_data.csv");

  ForestTrainer trainer = ForestTrainers::quantile_trainer(data, 10, quantiles);
  init_trainer(trainer);

  Forest forest = trainer.train(data);

  std::shared_ptr<PredictionStrategy> prediction_strategy(new QuantilePredictionStrategy(quantiles));
  ForestPredictor forest_predictor = create_forest_predictor(prediction_strategy);
  std::vector<std::vector<double>> predictions = forest_predictor.predict(forest, data);

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/quantile_test_predictions.csv");
  REQUIRE(predictions == expected_predictions);
}

TEST_CASE("causal forest predictions have not changed", "[causal, characterization]") {
  Data* data = loadDataFromFile("test/forest/resources/causal_test_data.csv");

  ForestTrainer trainer = ForestTrainers::instrumental_trainer(data, 10, 11, 11);
  init_trainer(trainer);

  Forest forest = trainer.train(data);

  std::shared_ptr<PredictionStrategy> prediction_strategy(new InstrumentalPredictionStrategy());
  ForestPredictor forest_predictor = create_forest_predictor(prediction_strategy);
  std::vector<std::vector<double>> predictions = forest_predictor.predict(forest, data);

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/causal_test_predictions.csv");

  REQUIRE(predictions.size() == expected_predictions.size());

  for (int i = 0; i < predictions.size(); i++) {
    std::vector<double> prediction = predictions[i];
    std::vector<double> expected_prediction = expected_predictions[i];

    REQUIRE(prediction.size() == 1);
    REQUIRE(expected_prediction.size() == 1);

    REQUIRE(equalDoubles(prediction[0], expected_prediction[0], 1e-2));
  }
}

TEST_CASE("regression forest predictions have not changed", "[regression, characterization]") {
  Data* data = loadDataFromFile("test/forest/resources/regression_test_data.csv");

  ForestTrainer trainer = ForestTrainers::regression_trainer(data, 10);
  init_trainer(trainer);

  Forest forest = trainer.train(data);

  std::shared_ptr<PredictionStrategy> prediction_strategy(new RegressionPredictionStrategy());
  ForestPredictor forest_predictor = create_forest_predictor(prediction_strategy);
  std::vector<std::vector<double>> predictions = forest_predictor.predict(forest, data);

  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::readCsvFile(
  "test/forest/resources/regression_test_predictions.csv");

  REQUIRE(predictions.size() == expected_predictions.size());

  for (int i = 0; i < predictions.size(); i++) {
    std::vector<double> prediction = predictions[i];
    std::vector<double> expected_prediction = expected_predictions[i];

    REQUIRE(prediction.size() == 1);
    REQUIRE(expected_prediction.size() == 1);

    REQUIRE(equalDoubles(prediction[0], expected_prediction[0], 1e-2));
  }
}