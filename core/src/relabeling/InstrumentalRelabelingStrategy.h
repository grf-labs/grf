#ifndef GRADIENTFOREST_INSTRUMENTALRELABELINGSTRATEGY_H
#define GRADIENTFOREST_INSTRUMENTALRELABELINGSTRATEGY_H

#include <unordered_map>
#include <vector>
#include "Observations.h"
#include "Tree.h"
#include "RelabelingStrategy.h"

class InstrumentalRelabelingStrategy: public RelabelingStrategy {
public:
  InstrumentalRelabelingStrategy();

  InstrumentalRelabelingStrategy(double split_regularization);

  std::unordered_map<size_t, double> relabel_outcomes(
      const Observations& observations,
      std::vector<size_t> &node_sampleIDs);

  DISALLOW_COPY_AND_ASSIGN(InstrumentalRelabelingStrategy);

private:
  double split_regularization;
};

#endif //GRADIENTFOREST_INSTRUMENTALRELABELINGSTRATEGY_H
