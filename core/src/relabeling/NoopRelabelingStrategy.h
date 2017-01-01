#ifndef GRADIENTFOREST_NOOPRELABELINGSTRATEGY_H
#define GRADIENTFOREST_NOOPRELABELINGSTRATEGY_H

#include "RelabelingStrategy.h"

class NoopRelabelingStrategy: public RelabelingStrategy {
public:
  std::unordered_map<size_t, double> relabel_outcomes(
      Observations *observations,
      std::vector<size_t> &node_sampleIDs);
};


#endif //GRADIENTFOREST_NOOPRELABELINGSTRATEGY_H
