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

  Author: Julie Tibshirani (jtibs@cs.stanford.edu)
 #-------------------------------------------------------------------------------*/

#include "Observations.h"

const std::string Observations::OUTCOME = "outcome";
const std::string Observations::TREATMENT = "treatment";
const std::string Observations::INSTRUMENT = "instrument";

Observations::Observations(const std::map<std::string, std::vector<double>>& observations_by_type,
                           size_t num_samples):
  observations_by_type(observations_by_type),
  num_samples(num_samples) {}

const std::vector<double>& Observations::get(std::string type) const {
  if (observations_by_type.find(type) == observations_by_type.end()) {
    throw std::runtime_error(
        "No observations of type " + type + " could be found.");
  }
  return observations_by_type.at(type);
}


