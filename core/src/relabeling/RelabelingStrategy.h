//
// Created by Julie Noelle Tibshirani on 10/30/16.
//
#ifndef GRADIENTFOREST_RELABELINGSTRATEGY_H
#define GRADIENTFOREST_RELABELINGSTRATEGY_H


#include <utility/Data.h>
#include <unordered_map>
#include <vector>

class RelabelingStrategy {
public:
  virtual std::unordered_map<size_t, double> relabelObservations(
      std::unordered_map<std::string, std::vector<double>> *observations,
      std::vector<size_t> &nodeSampleIDs) = 0;
};


#endif //GRADIENTFOREST_RELABELINGSTRATEGY_H
