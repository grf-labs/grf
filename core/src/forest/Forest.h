#ifndef GRADIENTFOREST_FOREST_H_
#define GRADIENTFOREST_FOREST_H_

#include "TreeTrainer.h"
#include "globals.h"
#include "Tree.h"
#include "Data.h"
#include "PredictionStrategy.h"
#include "Observations.h"

class Forest {
public:
  static Forest create(std::vector<std::shared_ptr<Tree>> trees,
                       Data* data,
                       std::unordered_map<std::string, size_t> observables);
  Forest(std::vector<std::shared_ptr<Tree>> trees,
         Observations observations);

  const Observations get_observations() const {
    return observations;
  };

  const std::vector<std::shared_ptr<Tree>> get_trees() const {
    return trees;
  }

  const std::vector<std::vector<size_t>> get_inbag_counts() const {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : trees) {
      result.push_back(tree->get_inbag_counts());
    }
    return result;
  }

protected:
  std::vector<std::shared_ptr<Tree>> trees;
  Observations observations;
};

#endif /* GRADIENTFOREST_FOREST_H_ */
