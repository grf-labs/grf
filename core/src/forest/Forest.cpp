#include "Data.h"
#include "Forest.h"

Forest Forest::create(std::vector<std::shared_ptr<Tree>> trees,
                      Data* data,
                      std::unordered_map<std::string, size_t> observables) {
  std::map<std::string, std::vector<double>> observations_by_type;
  size_t num_samples = data->getNumRows();
  for (auto it : observables) {
    std::string type = it.first;
    size_t index = it.second;
    for (int row = 0; row < num_samples; row++) {
      observations_by_type[type].push_back(data->get(row, index));
    }
  }

  Observations observations(observations_by_type, num_samples);
  return Forest(trees, observations);
}

Forest::Forest(const std::vector<std::shared_ptr<Tree>>& trees,
               const Observations& observations):
  trees(trees), observations(observations) {}
