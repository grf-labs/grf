#include "Observations.h"

const std::string Observations::OUTCOME = "outcome";
const std::string Observations::TREATMENT = "treatment";
const std::string Observations::INSTRUMENT = "instrument";

Observations::Observations(std::map<std::string, std::vector<double>> observations_by_type,
                           size_t num_samples):
  observations_by_type(observations_by_type),
  num_samples(num_samples) {}

const std::vector<double> Observations::get(std::string type) const {
  if (observations_by_type.find(type) == observations_by_type.end()) {
    throw std::runtime_error(
        "No observations of type " + type + " could be found.");
  }
  return observations_by_type.at(type);
}


