#ifndef GRADIENTFOREST_FORESTMODEL_H
#define GRADIENTFOREST_FORESTMODEL_H

#include "RelabelingStrategy.h"
#include "SplittingRule.h"
#include "PredictionStrategy.h"

#include "Tree.h"
#include "TreeTrainer.h"
#include "Forest.h"

#include "ForestTrainer.h"
#include "ForestPredictor.h"

#include <thread>
#include <future>

class ForestModel {
public:
  ForestModel(ForestTrainer &trainer, ForestPredictor &predictor);

  Forest* train(Data* data);
  std::vector<std::vector<double>> predict(Forest* forest, Data* prediction_data);

private:
  ForestTrainer trainer;
  ForestPredictor predictor;
};


#endif //GRADIENTFOREST_FORESTMODEL_H
