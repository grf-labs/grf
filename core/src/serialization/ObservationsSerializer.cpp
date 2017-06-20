/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "serialization/ObservationsSerializer.h"
#include "commons/utility.h"

void ObservationsSerializer::serialize(std::ostream& stream, const Observations& observations) {
  size_t num_samples = observations.get_num_samples();
  stream.write((char*) &num_samples, sizeof(num_samples));

  const auto& observations_by_type = observations.get_observations_by_type();
  size_t num_types = observations_by_type.size();
  stream.write((char*) &num_types, sizeof(num_types));

  for (auto& observations : observations_by_type) {
    write_vector(observations, stream);
  }
}

Observations ObservationsSerializer::deserialize(std::istream& stream) {
  size_t num_samples;
  stream.read((char*) &num_samples, sizeof(num_samples));

  size_t num_types;
  stream.read((char*) &num_types, sizeof(num_types));

  std::vector<std::vector<double>> observations_by_type(num_types);
  for (size_t i = 0; i < num_types; ++i) {
    read_vector(observations_by_type[i], stream);
  }

  return Observations(observations_by_type, num_samples);
}
