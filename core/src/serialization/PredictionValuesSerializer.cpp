#include "PredictionValuesSerializer.h"
#include "utility.h"

void PredictionValuesSerializer::serialize(std::ostream& stream, const PredictionValues& values) {
  size_t num_nodes = values.get_num_nodes();
  stream.write((char*) &num_nodes, sizeof(num_nodes));

  auto values_by_type = values.get_values_by_type();
  size_t num_types = values_by_type.size();
  stream.write((char*) &num_types, sizeof(num_types));

  for (auto it = values_by_type.begin(); it != values_by_type.end(); it++) {
    std::string type = it->first;
    stream.write((char*) &type, sizeof(type));
    saveVector1D(it->second, stream);
  }
}

PredictionValues PredictionValuesSerializer::deserialize(std::istream& stream) {
  size_t num_nodes;
  stream.read((char*) &num_nodes, sizeof(num_nodes));

  size_t num_types;
  stream.read((char*) &num_types, sizeof(num_types));

  std::map<std::string, std::vector<double>> values_by_type;
  for (int i = 0; i < num_types; i++) {
    std::string type;
    stream.read((char*) &type, sizeof(type));

    std::vector<double> values;
    readVector1D(values, stream);
    values_by_type[type] = values;
  }

  return PredictionValues(values_by_type, num_nodes);
}
