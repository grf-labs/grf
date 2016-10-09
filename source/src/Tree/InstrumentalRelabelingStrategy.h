#ifndef RANGER_INSTRUMENTALRELABELINGSTRATEGY_H
#define RANGER_INSTRUMENTALRELABELINGSTRATEGY_H

#include <unordered_map>
#include "Tree.h"

class InstrumentalRelabelingStrategy {
public:

  InstrumentalRelabelingStrategy(Data* data,
                                 size_t dependent_varID,
                                 size_t treatment_varID,
                                 size_t instrument_varID);
  std::unordered_map<size_t, double> relabelResponses(std::vector<size_t>& nodeSampleIDs);

private:
  bool equalDoubles(double first, double second);
  int sgn(double val);

  Data* data;
  size_t dependent_varID;
  size_t treatment_varID;
  size_t instrument_varID;

  DISALLOW_COPY_AND_ASSIGN(InstrumentalRelabelingStrategy);

};

#endif //RANGER_INSTRUMENTALRELABELINGSTRATEGY_H
