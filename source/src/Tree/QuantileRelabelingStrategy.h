//
// Created by Julie Noelle Tibshirani on 10/9/16.
//

#ifndef RANGER_QUANTILERELABELINGSTRATEGY_H
#define RANGER_QUANTILERELABELINGSTRATEGY_H

#include "TreeFactory.h"

class QuantileRelabelingStrategy {
public:
  QuantileRelabelingStrategy(std::vector<double>* quantiles);
  std::vector<uint>* relabelResponses(std::vector<double> *responses);

private:
  std::vector<double>* quantiles;
};


#endif //RANGER_QUANTILERELABELINGSTRATEGY_H
