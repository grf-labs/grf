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
  static Forest* create(std::vector<Tree*>* trees,
                        Data* data,
                        std::unordered_map<std::string, size_t> observables);
  Forest(std::vector<Tree*>* trees,
         Observations observations);

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
      result.push_back(tree->get_inbag_counts());
    }
    return result;
  }

protected:
  std::vector<Tree*>* trees;
  Observations observations;

private:
  DISALLOW_COPY_AND_ASSIGN(Forest);
};

#endif /* GRADIENTFOREST_FOREST_H_ */
