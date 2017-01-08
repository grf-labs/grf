#ifndef GRADIENTFOREST_TREESERIALIZER_H
#define GRADIENTFOREST_TREESERIALIZER_H


class TreeSerializer {
public:
  void serialize(std::ostream& stream, std::shared_ptr<Tree> tree);
  std::shared_ptr<Tree> deserialize(std::istream& stream);
};


#endif //GRADIENTFOREST_TREESERIALIZER_H
