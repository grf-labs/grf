//
// Created by Julie Noelle Tibshirani on 10/9/16.
//

#ifndef RANGER_QUANTILERELABELINGSTRATEGY_H
#define RANGER_QUANTILERELABELINGSTRATEGY_H

#include "Tree.h"
#include "RelabelingStrategy.h"

class QuantileRelabelingStrategy: public RelabelingStrategy {
public:
  QuantileRelabelingStrategy(std::vector<double>* quantiles);
  std::unordered_map<size_t, double> relabelObservations(
      std::unordered_map<std::string, std::vector<double>> *observations,
      std::vector<size_t> &node_sampleIDs);
private:
  std::vector<double>* quantiles;
};


#endif //RANGER_QUANTILERELABELINGSTRATEGY_H
