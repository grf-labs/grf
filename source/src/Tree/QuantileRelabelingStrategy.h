//
// Created by Julie Noelle Tibshirani on 10/9/16.
//

#ifndef RANGER_QUANTILERELABELINGSTRATEGY_H
#define RANGER_QUANTILERELABELINGSTRATEGY_H

#include "TreeFactory.h"

class QuantileRelabelingStrategy {
public:
  QuantileRelabelingStrategy(std::vector<double>* quantiles, size_t dependent_varID);
  std::vector<uint>* relabelResponses(Data* data,
                                      std::vector<size_t>& nodeSampleIDs);

private:
  std::vector<double>* quantiles;
  size_t dependent_varID;
};


#endif //RANGER_QUANTILERELABELINGSTRATEGY_H
