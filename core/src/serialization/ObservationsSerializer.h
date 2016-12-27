#ifndef GRADIENTFOREST_OBSERVATIONSSERIALIZER_H
#define GRADIENTFOREST_OBSERVATIONSSERIALIZER_H


#include <ostream>
#include <istream>
#include "Observations.h"

class ObservationsSerializer {
public:
  void serialize(std::ostream& stream, Observations* observations);
  Observations* deserialize(std::istream& stream);
};


#endif //GRADIENTFOREST_OBSERVATIONSSERIALIZER_H
