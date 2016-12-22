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
  std::unordered_map<size_t, double> relabel_outcomes(
      Observations *observations,
      std::vector<size_t> &node_sampleIDs);

  DISALLOW_COPY_AND_ASSIGN(InstrumentalRelabelingStrategy);

};

#endif //GRADIENTFOREST_INSTRUMENTALRELABELINGSTRATEGY_H
