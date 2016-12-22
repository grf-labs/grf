#ifndef GRADIENTFOREST_QUANTILERELABELINGSTRATEGY_H
#define GRADIENTFOREST_QUANTILERELABELINGSTRATEGY_H

#include "Observations.h"
#include "Tree.h"
#include "RelabelingStrategy.h"

class QuantileRelabelingStrategy: public RelabelingStrategy {
public:
  QuantileRelabelingStrategy(std::vector<double>* quantiles);
  std::unordered_map<size_t, double> relabel_outcomes(
      Observations *observations,
      std::vector<size_t> &node_sampleIDs);
private:
  std::vector<double>* quantiles;
};


#endif //GRADIENTFOREST_QUANTILERELABELINGSTRATEGY_H
