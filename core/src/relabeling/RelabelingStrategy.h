#ifndef GRADIENTFOREST_RELABELINGSTRATEGY_H
#define GRADIENTFOREST_RELABELINGSTRATEGY_H


#include "Data.h"
#include "Observations.h"
#include <unordered_map>
#include <vector>

class RelabelingStrategy {
public:
  virtual std::unordered_map<size_t, double> relabel_outcomes(
      const Observations& observations,
      const std::vector<size_t>& nodeSampleIDs) = 0;
};


#endif //GRADIENTFOREST_RELABELINGSTRATEGY_H
