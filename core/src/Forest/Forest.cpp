
#include <utility/Data.h>
#include "Forest.h"

Forest::Forest(std::vector<Tree *>* trees,
               Data *data,
               std::unordered_map<std::string, size_t> observables) :
    trees(trees),
    data(data),
    observables_by_ID(observables) {
  for (auto it : observables) {
    std::string name = it.first;
    size_t index = it.second;

    for (int row = 0; row < data->getNumRows(); row++) {
      original_observations[name].push_back(data->get(row, index));
    }
  }
}

Forest::~Forest() {
  for (auto& tree : *trees) {
    delete tree;
  }
}
