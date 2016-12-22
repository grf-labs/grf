#ifndef GRADIENTFOREST_FOREST_H_
#define GRADIENTFOREST_FOREST_H_

#include "TreeModel.h"
#include "globals.h"
#include "Tree.h"
#include "Data.h"
#include "PredictionStrategy.h"
#include "Observations.h"

class Forest {
public:
  Forest(std::vector<Tree*>* trees,
         Data* data,
         std::unordered_map<std::string, size_t> observables);
  virtual ~Forest();

  const Observations get_observations() const {
    return observations;
  };

  const std::vector<Tree*>* get_trees() const {
    return trees;
  }

  const std::vector<std::vector<size_t>> get_inbag_counts() const {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : *trees) {
      result.push_back(tree->getInbagCounts());
    }
    return result;
  }

protected:
  std::vector<Tree*>* trees;

  Data* data;
  Observations observations;

private:
  DISALLOW_COPY_AND_ASSIGN(Forest);
};

#endif /* GRADIENTFOREST_FOREST_H_ */
