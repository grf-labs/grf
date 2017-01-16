#ifndef GRADIENTFOREST_PREDICTIONVALUESSERIALIZER_H
#define GRADIENTFOREST_PREDICTIONVALUESSERIALIZER_H

#include "PredictionValues.h"

class PredictionValuesSerializer {
public:
  void serialize(std::ostream& stream, const PredictionValues& observations);
  PredictionValues deserialize(std::istream& stream);
};


#endif //GRADIENTFOREST_PREDICTIONVALUESSERIALIZER_H
