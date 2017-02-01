#include "PredictionStrategy.h"
#include "utility.h"
#include "ForestPredictor.h"
#include "ForestTrainer.h"
#include "FileTestUtilities.h"
#include "ForestTestUtilities.h"

#include "ForestTrainers.h"
#include "ForestPredictors.h"

#include "catch.hpp"


TEST_CASE("honest regression forests are shift invariant", "[regression, forest]") {
  // Run the original forest.
  Data* data = loadDataFromFile("test/forest/resources/gaussian_data.csv");
  uint outcome_index = 10;

  ForestTrainer trainer = ForestTrainers::regression_trainer(data, outcome_index);
  ForestTestUtilities::init_honest_trainer(trainer);

  Forest forest = trainer.train(data);
  ForestPredictor predictor = ForestPredictors::regression_predictor(4);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  // Shift each outcome by 1, and re-run the forest.
  bool error;
  for (size_t r = 0; r < data->getNumRows(); r++) {
    double outcome = data->get(r, outcome_index);
    data->set(outcome_index, r, outcome + 1, error);
  }

  ForestTrainer shifted_trainer = ForestTrainers::regression_trainer(data, outcome_index);
  ForestTestUtilities::init_trainer(shifted_trainer);

  Forest shifted_forest = trainer.train(data);
  ForestPredictor shifted_predictor = ForestPredictors::regression_predictor(4);
  std::vector<Prediction> shifted_predictions = shifted_predictor.predict_oob(shifted_forest, data);

  REQUIRE(predictions.size() == shifted_predictions.size());
  double delta = 0.0;
  for (size_t i = 0; i < predictions.size(); i++) {
    Prediction prediction = predictions[i];
    Prediction shifted_prediction = shifted_predictions[i];

    double value = prediction.get_predictions()[0];
    double shifted_value = shifted_prediction.get_predictions()[0];

    delta += shifted_value - value;
  }

  REQUIRE(equalDoubles(delta / predictions.size(), 1, 1e-1));
  
  delete data;
}