#ifndef RANGER_INSTRUMENTALRELABELINGSTRATEGY_H
#define RANGER_INSTRUMENTALRELABELINGSTRATEGY_H

#include <unordered_map>
#include "TreeFactory.h"
#include "RelabelingStrategy.h"

class InstrumentalRelabelingStrategy: public RelabelingStrategy {
public:
  InstrumentalRelabelingStrategy(std::unordered_map<std::string, size_t> observables);
  std::unordered_map<size_t, double> relabelResponses(Data* data,
                                                      std::vector<size_t>& nodeSampleIDs);

private:
  bool equalDoubles(double first, double second);
  int sgn(double val);

  size_t dependent_varID;
  size_t treatment_varID;
  size_t instrument_varID;

  DISALLOW_COPY_AND_ASSIGN(InstrumentalRelabelingStrategy);

};

#endif //RANGER_INSTRUMENTALRELABELINGSTRATEGY_H
