#ifndef GRADIENTFOREST_TREEPREDICTOR_H
#define GRADIENTFOREST_TREEPREDICTOR_H

#include "Data.h"
#include "Tree.h"

class TreePredictor {
public:
  std::vector<size_t> get_terminal_nodeIDs(std::shared_ptr<Tree>,
                                           const Data* prediction_data,
                                           std::vector<size_t> sampleIDs);
};


#endif //GRADIENTFOREST_TREEPREDICTOR_H
