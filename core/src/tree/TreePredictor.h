#ifndef GRADIENTFOREST_TREEPREDICTOR_H
#define GRADIENTFOREST_TREEPREDICTOR_H

#include <memory>

#include "Data.h"
#include "Tree.h"

class TreePredictor {
public:
  std::vector<size_t> get_terminal_nodeIDs(std::shared_ptr<Tree>,
                                           Data* prediction_data,
                                           const std::vector<size_t>& sampleIDs);
};


#endif //GRADIENTFOREST_TREEPREDICTOR_H
