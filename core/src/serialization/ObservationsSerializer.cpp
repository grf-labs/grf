/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

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
    write_string(type, stream);
    write_vector(it->second, stream);
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
    read_string(type, stream);
    read_vector(observations_by_type[type], stream);
  }

  return Observations(observations_by_type, num_samples);
}
