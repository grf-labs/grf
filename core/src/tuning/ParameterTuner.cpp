/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include <limits.h>
#include "ParameterTuner.h"

ParameterTuner::ParameterTuner(const ForestTrainer& trainer,
                               const ForestPredictor& predictor,
                               uint outcome_index):
    trainer(trainer), predictor(predictor), outcome_index(outcome_index) {}

uint ParameterTuner::tune_min_node_size(Data* data,
                                        ForestOptions& options) {
  uint upper_bound = (uint) (data->get_num_rows() * options.get_sample_fraction()) / 4;

  uint min_node_size = 10;
  double best_mse = std::numeric_limits<double>::max();
  uint best_min_node_size = NAN;

  while (min_node_size < upper_bound) {
    options.set_min_node_size(min_node_size);

    const Forest forest = trainer.train(data, options);
    std::vector<Prediction> predictions = predictor.predict_oob(forest, data);
    double mse = calculate_mse(predictions, data);

    if (mse < best_mse) {
      best_mse = mse;
      best_min_node_size = min_node_size;
    }
    min_node_size *= 2;
  }

  return best_min_node_size;
}

double ParameterTuner::calculate_mse(const std::vector<Prediction>& predictions,
                                     Data* data) {
  double mean_squared_error = 0;
  for (size_t i = 0; i < predictions.size(); ++i) {
    double actual_outcome = data->get(i, outcome_index);
    const Prediction& prediction = predictions[i];
    double predicted_outcome = prediction.get_predictions().at(0);

    double difference = (actual_outcome - predicted_outcome);
    mean_squared_error += difference * difference;
  }
  return mean_squared_error / predictions.size();
}
