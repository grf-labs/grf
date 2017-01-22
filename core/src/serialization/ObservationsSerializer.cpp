#include "ObservationsSerializer.h"
#include "utility.h"

void ObservationsSerializer::serialize(std::ostream& stream, const Observations& observations) {
  size_t num_samples = observations.get_num_samples();
  stream.write((char*) &num_samples, sizeof(num_samples));

  const auto& observations_by_type = observations.get_observations_by_type();
  size_t num_types = observations_by_type.size();
  stream.write((char*) &num_types, sizeof(num_types));

  for (auto it = observations_by_type.begin(); it != observations_by_type.end(); it++) {
    std::string type = it->first;
    saveString(type, stream);
    saveVector1D(it->second, stream);
  }
}

Observations ObservationsSerializer::deserialize(std::istream& stream) {
  size_t num_samples;
  stream.read((char*) &num_samples, sizeof(num_samples));

  size_t num_types;
  stream.read((char*) &num_types, sizeof(num_types));

  std::map<std::string, std::vector<double>> observations_by_type;
  for (size_t i = 0; i < num_types; i++) {
    std::string type;
    readString(type, stream);
    readVector1D(observations_by_type[type], stream);
  }

  return Observations(observations_by_type, num_samples);
}
