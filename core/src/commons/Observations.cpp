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

#include <unordered_map>
#include "Observations.h"

const std::size_t Observations::OUTCOME = 0;
const std::size_t Observations::TREATMENT = 1;
const std::size_t Observations::INSTRUMENT = 2;

Observations::Observations():
  observations_by_type(std::vector<std::vector<double>>()),
  num_samples(0) {}

Observations::Observations(const std::vector<std::vector<double>>& observations_by_type,
                           size_t num_samples):
  observations_by_type(observations_by_type),
  num_samples(num_samples) {}

double Observations::get(std::size_t type, size_t sample) const {
  return observations_by_type.at(type).at(sample);
}

const std::vector<std::vector<double>>& Observations::get_observations_by_type() const {
  return observations_by_type;
}

size_t Observations::get_num_samples() const {
  return num_samples;
}

size_t Observations::get_num_types() const {
  return observations_by_type.size();
}

