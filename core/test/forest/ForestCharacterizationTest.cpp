#include <MacTypes.h>
#include "PredictionStrategy.h"
#include "utility.h"
#include "ForestPredictor.h"
#include "ForestTrainer.h"
#include "FileTestUtilities.h"

#include "ForestTrainers.h"
#include "ForestPredictors.h"

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
  double sample_fraction = 0.7;
  bool honesty = false;
  uint ci_bag_size = 2;

  trainer.init(mtry, num_trees, verbose_out, seed, num_threads, load_forest_filename,
               min_node_size, no_split_variables, split_select_weights_file, always_split_variable_names,
               sample_with_replacement, memory_saving_splitting, case_weights_file, sample_fraction,
               honesty, ci_bag_size);
}

bool equal_predictions(std::vector<std::vector<double>> predictions,
                       std::vector<std::vector<double>> expected_predictions) {
  if (predictions.size() != expected_predictions.size()) {
    return false;
  }
  for (int i = 0; i < predictions.size(); i++) {
    std::vector<double> prediction = predictions[i];
    std::vector<double> expected_prediction = expected_predictions[i];
    if (prediction.size() != expected_prediction.size()) {
      return false;
    }

    for (int j = 0; j < prediction.size(); j++) {
      if (!equalDoubles(prediction[j], expected_prediction[j], 1e-2)) {
        return false;
      }
    }
  }
  return true;
}

TEST_CASE("quantile forest predictions have not changed", "[quantile, characterization]") {
  std::vector<double> quantiles({0.25, 0.5, 0.75});
  Data *data = loadDataFromFile("test/forest/resources/quantile_data.csv");

  ForestTrainer trainer = ForestTrainers::quantile_trainer(data, 10, quantiles);
  init_trainer(trainer);
  Forest forest = trainer.train(data);

  ForestPredictor predictor = ForestPredictors::quantile_predictor(4, quantiles);

  std::vector<std::vector<double>> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/quantile_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> predictions = predictor.predict(forest, data);
  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/quantile_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("causal forest predictions have not changed", "[causal, characterization]") {
  Data* data = loadDataFromFile("test/forest/resources/causal_data.csv");

  ForestTrainer trainer = ForestTrainers::instrumental_trainer(data, 10, 11, 11, 0);
  init_trainer(trainer);

  Forest forest = trainer.train(data);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(4);

  std::vector<std::vector<double>> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/causal_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> predictions = predictor.predict(forest, data);
  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/causal_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}

TEST_CASE("regression forest predictions have not changed", "[regression, characterization]") {
  Data* data = loadDataFromFile("test/forest/resources/regression_data.csv");

  ForestTrainer trainer = ForestTrainers::regression_trainer(data, 10);
  init_trainer(trainer);

  Forest forest = trainer.train(data);

  ForestPredictor predictor = ForestPredictors::regression_predictor(4);

  std::vector<std::vector<double>> oob_predictions = predictor.predict_oob(forest, data);
  std::vector<std::vector<double>> expected_oob_predictions = FileTestUtilities::readCsvFile(
      "test/forest/resources/regression_oob_predictions.csv");
  REQUIRE(equal_predictions(oob_predictions, expected_oob_predictions));

  std::vector<std::vector<double>> predictions = predictor.predict(forest, data);
  std::vector<std::vector<double>> expected_predictions = FileTestUtilities::readCsvFile(
  "test/forest/resources/regression_predictions.csv");
  REQUIRE(equal_predictions(predictions, expected_predictions));
}
