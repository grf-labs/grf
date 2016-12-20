#include <prediction/PredictionStrategy.h>
#include <utility/utility.h>
#include <forest/ForestModel.h>
#include <relabeling/QuantileRelabelingStrategy.h>
#include <prediction/QuantilePredictionStrategy.h>
#include <splitting/RegressionSplittingRule.h>
#include <splitting/ProbabilitySplittingRule.h>
#include <relabeling/InstrumentalRelabelingStrategy.h>
#include <prediction/InstrumentalPredictionStrategy.h>

#include "catch.hpp"

// TODO(jtibs): Clean up these helper methods.
void initializeForestModel(ForestModel* forest_model) {
  uint mtry = 3;
  uint num_trees = 4;
  std::ostream* verbose_out = &std::cout;
  uint seed = 42;
  uint num_threads = 1;
  std::string load_forest_filename = "";
  uint min_node_size = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  std::string split_select_weights_file = "";
  std::vector<std::string> always_split_variable_names;
  bool sample_with_replacement = true;
  bool memory_saving_splitting = false;
  std::string case_weights_file = "";
  double sample_fraction = 1;

  forest_model->initCpp(mtry, num_trees, verbose_out, seed, num_threads, load_forest_filename,
                        min_node_size, split_select_weights_file, always_split_variable_names, sample_with_replacement,
                        memory_saving_splitting, case_weights_file, sample_fraction);
}

std::vector<std::vector<double>> loadPredictionsFile(std::string file_name) {
  std::ifstream predictions_file = std::ifstream();
  predictions_file.open(file_name, std::ios::binary);

  std::string delimiter = " ";
  std::vector<std::vector<double>> predictions;
  std::string line;

  while (std::getline(predictions_file, line)) {
    size_t position = 0;
    std::string token;

    std::vector<double> prediction;
    while ((position = line.find(delimiter)) != std::string::npos) {
      token = line.substr(0, position);
      prediction.push_back(std::stod(token));
      line.erase(0, position + delimiter.length());
    }

    if (position != line.length()) {
      token = line.substr(0, line.length());
      prediction.push_back(std::stod(token));
    }

    predictions.push_back(prediction);
  }

  predictions_file.close();
  return predictions;
}

TEST_CASE("quantile forest predictions have not changed", "[quantile, characterization]") {
  std::vector<double>* quantiles = new std::vector<double>({0.25, 0.5, 0.75});
  std::unordered_map<std::string, size_t> observables = {{"outcome", 10}};
  Data* data = loadDataFromFile("test/forest/resources/quantile_test_data.csv");

  RelabelingStrategy* relabeling_strategy = new QuantileRelabelingStrategy(quantiles);
  SplittingRule* splitting_rule = new ProbabilitySplittingRule(data, quantiles->size());
  PredictionStrategy* prediction_strategy = new QuantilePredictionStrategy(quantiles);

  ForestModel* forest_model = new ForestModel(observables,
      relabeling_strategy,
      splitting_rule,
      prediction_strategy);
  initializeForestModel(forest_model);

  Forest* forest = forest_model->train(data);
  std::vector<std::vector<double>> predictions = forest_model->predict(forest, data);

  std::vector<std::vector<double>> expected_predictions = loadPredictionsFile(
      "test/forest/resources/quantile_test_predictions.csv");
  REQUIRE(predictions == expected_predictions);
}

TEST_CASE("causal forest predictions have not changed", "[causal, characterization]") {
  std::unordered_map<std::string, size_t> observables = {
      {"outcome", 10}, {"treatment", 11}, {"instrument", 11}};
  Data* data = loadDataFromFile("test/forest/resources/causal_test_data.csv");

  RelabelingStrategy* relabeling_strategy = new InstrumentalRelabelingStrategy();
  SplittingRule* splitting_rule = new RegressionSplittingRule(data);
  PredictionStrategy* prediction_strategy = new InstrumentalPredictionStrategy();

  ForestModel* forest_model = new ForestModel(observables,
                                              relabeling_strategy,
                                              splitting_rule,
                                              prediction_strategy);
  initializeForestModel(forest_model);

  Forest* forest = forest_model->train(data);
  std::vector<std::vector<double>> predictions = forest_model->predict(forest, data);

  std::vector<std::vector<double>> expected_predictions = loadPredictionsFile(
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