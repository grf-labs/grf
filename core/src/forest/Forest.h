#ifndef GRADIENTFOREST_FOREST_H_
#define GRADIENTFOREST_FOREST_H_

#include <memory>

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
  Forest(const std::vector<std::shared_ptr<Tree>>& trees,
         const Observations& observations);

  const Observations& get_observations() const {
    return observations;
  };

  const std::vector<std::shared_ptr<Tree>>& get_trees() const {
    return trees;
  }

protected:
  std::vector<std::shared_ptr<Tree>> trees;
  Observations observations;
};

#endif /* GRADIENTFOREST_FOREST_H_ */
