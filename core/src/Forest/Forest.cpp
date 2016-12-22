
#include "Data.h"
#include "Forest.h"

Forest::Forest(std::vector<Tree *>* trees,
               Data *data,
               std::unordered_map<std::string, size_t> observables) :
    trees(trees),
    data(data) {

  std::unordered_map<std::string, std::vector<double>> observationsByType;
  size_t num_samples = data->getNumRows();
  for (auto it : observables) {
    std::string type = it.first;
    size_t index = it.second;
    for (int row = 0; row < num_samples; row++) {
      observationsByType[type].push_back(data->get(row, index));
    }
  }
  this->observations = Observations(observationsByType, num_samples);
}

Forest::~Forest() {
  for (auto& tree : *trees) {
    delete tree;
  }
}
