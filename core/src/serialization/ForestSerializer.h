#ifndef GRADIENTFOREST_FORESTSERIALIZER_H
#define GRADIENTFOREST_FORESTSERIALIZER_H


#include "Forest.h"
#include "ObservationsSerializer.h"
#include "TreeSerializer.h"

class ForestSerializer {
public:
  void serialize(std::ostream& stream, Forest* forest);
  Forest* deserialize(std::istream& stream);

private:
  ObservationsSerializer observation_serializer;
  TreeSerializer tree_serializer;
};


#endif //GRADIENTFOREST_FORESTSERIALIZER_H
