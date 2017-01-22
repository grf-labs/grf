#include "PredictionValues.h"

const std::string PredictionValues::AVERAGE = "average";

PredictionValues::PredictionValues() {}

PredictionValues::PredictionValues(const std::map<std::string, std::vector<double>>& values_by_type,
                                   size_t num_nodes):
  values_by_type(values_by_type),
  num_nodes(num_nodes) {}

const std::vector<double>& PredictionValues::get(std::string type) const {
  if (values_by_type.find(type) == values_by_type.end()) {
    throw std::runtime_error(
        "No observations of type " + type + " could be found.");
  }
  return values_by_type.at(type);
}


