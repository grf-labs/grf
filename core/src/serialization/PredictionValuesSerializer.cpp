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

#include "PredictionValuesSerializer.h"
#include "utility.h"

void PredictionValuesSerializer::serialize(std::ostream& stream, const PredictionValues& prediction_values) {
  size_t num_nodes = prediction_values.get_num_nodes();
  stream.write((char*) &num_nodes, sizeof(num_nodes));

  auto values_by_type = prediction_values.get_values_by_type();
  size_t num_types = values_by_type.size();
  stream.write((char*) &num_types, sizeof(num_types));

  for (auto& values : values_by_type) {
    write_vector(values, stream);
  }
}

PredictionValues PredictionValuesSerializer::deserialize(std::istream& stream) {
  size_t num_nodes;
  stream.read((char*) &num_nodes, sizeof(num_nodes));

  size_t num_types;
  stream.read((char*) &num_types, sizeof(num_types));

  std::vector<std::vector<double>> values_by_type(num_types);
  for (size_t i = 0; i < num_types; i++) {
    read_vector(values_by_type[i], stream);
  }

  return PredictionValues(values_by_type, num_nodes);
}
