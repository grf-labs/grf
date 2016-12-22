#include "Observations.h"

Observations::Observations(std::unordered_map<std::string, std::vector<double>> observationsByType,
                           size_t size):
  observationsByType(observationsByType),
  num_samples(size) {}

std::vector<double> Observations::get(std::string type) {
  if (observationsByType.find(type) == NULL) {
    throw std::runtime_error(
        "No observations of type " + type + " could be found.");
  }
  return observationsByType[type];
}


