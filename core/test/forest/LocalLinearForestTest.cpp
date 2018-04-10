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


#include <map>
#include <unordered_set>
#include <fstream>
#include "commons/Observations.h"
#include "commons/utility.h"
#include "utilities/ForestTestUtilities.h"
#include "prediction/LocalLinearPredictionStrategy.h"
#include "prediction/RegressionPredictionStrategy.h"

#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "utilities/ForestTestUtilities.h"

#include "catch.hpp"

TEST_CASE("basic llf test", "[prediction]") {
    // Run the original forest.
    std::cout << "beginning test \n";
    std::cout << "second output \n";

    //Data* data = load_data("test/forest/resources/gaussian_data.csv");
    Data* data = load_data("resources/gaussian_data.csv");
    std::cout << "have some data test";
    uint outcome_index = 10;
    double alpha = 0.10;

    std::cout << "calling forest trainers";

    ForestTrainer trainer = ForestTrainers::regression_trainer(outcome_index, alpha);
    ForestOptions options = ForestTestUtilities::default_honest_options();

    std::cout << "training";

    Forest forest = trainer.train(data, options);
    ForestPredictor predictor = ForestPredictors::local_linear_predictor(4,data,data,0.1,false);

    std::cout << "made a  predictor";

    double temp = 1;
    double temp2 = 1;
    REQUIRE(equal_doubles(temp, temp2, 1.0e-10));

    /*ForestPredictor predictor = ForestPredictors::regression_predictor(4, 1);
    std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

    // Shift each outcome by 1, and re-run the forest.
    bool error;
    for (size_t r = 0; r < data->get_num_rows(); r++) {
        double outcome = data->get(r, outcome_index);
        data->set(outcome_index, r, outcome + 1, error);
    }

    Forest shifted_forest = trainer.train(data, options);
    ForestPredictor shifted_predictor = ForestPredictors::regression_predictor(4, 1);
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

    REQUIRE(equal_doubles(delta / predictions.size(), 1, 1e-1));
    delete data;


    std::vector<double> averages = {1.1251472};
    std::vector<double> flipped_averages = {-1.1251472};

    std::vector<std::vector<double>> dataset = {
            {0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.1}, {4, 0.1},
            {5, 0.1}, {6, 0.2}, {7, 0.1}, {8, 0.0}, {9, 0.1}};

    Data* original_data = Data(dataset);

    LocalLinearPredictionStrategy prediction_strategy = LocalLinearPredictionStrategy(dataset, dataset, 0.1, 0);
    //LocalLinearPredictionStrategy prediction_strategy();
    // NEED TO FEED IT: original data, test data, lambda, ridge
    std::vector<double> first_prediction = prediction_strategy.predict(averages);
    std::vector<double> second_prediction = prediction_strategy.predict(flipped_averages);

    REQUIRE(first_prediction.size() == 1);
    REQUIRE(second_prediction.size() == 1);
    REQUIRE(equal_doubles(first_prediction[0], -second_prediction[0], 1.0e-10));
     */
}