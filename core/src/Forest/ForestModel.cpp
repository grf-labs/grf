#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <ctime>
#include <thread>
#include <future>
#include "utility.h"
#include "ForestModel.h"

ForestModel::ForestModel(ForestTrainer &trainer, ForestPredictor &predictor):
  trainer(trainer), predictor(predictor) {}

Forest* ForestModel::train(Data* data) {
  return trainer.train(data);
}

std::vector<std::vector<double>> ForestModel::predict(Forest* forest, Data* prediction_data) {
  return predictor.predict(forest, prediction_data);
}