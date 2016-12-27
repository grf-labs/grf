#ifndef GRADIENTFOREST_TREESERIALIZER_H
#define GRADIENTFOREST_TREESERIALIZER_H


class TreeSerializer {
public:
  void serialize(std::ostream& stream, Tree* tree);
  Tree* deserialize(std::istream& stream);
};


#endif //GRADIENTFOREST_TREESERIALIZER_H
