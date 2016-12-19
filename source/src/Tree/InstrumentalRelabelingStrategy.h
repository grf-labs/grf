#ifndef RANGER_INSTRUMENTALRELABELINGSTRATEGY_H
#define RANGER_INSTRUMENTALRELABELINGSTRATEGY_H

#include <unordered_map>
#include <vector>
#include "Tree.h"
#include "RelabelingStrategy.h"

class InstrumentalRelabelingStrategy: public RelabelingStrategy {
public:
  std::unordered_map<size_t, double> relabelObservations(
      std::unordered_map<std::string, std::vector<double>> *observations,
      std::vector<size_t> &node_sampleIDs);

private:
  bool equalDoubles(double first, double second);
  int sgn(double val);

  DISALLOW_COPY_AND_ASSIGN(InstrumentalRelabelingStrategy);

};

#endif //RANGER_INSTRUMENTALRELABELINGSTRATEGY_H
