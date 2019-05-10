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

#include "prediction/PredictionValues.h"

PredictionValues::PredictionValues():
  num_nodes(0),
  num_types(0) {}

PredictionValues::PredictionValues(const std::vector<std::vector<double>>& values,
                                   size_t num_types):
  values(values),
  num_nodes(values.size()),
  num_types(num_types) {}



double PredictionValues::get(std::size_t node, size_t type) const {
  return values.at(node).at(type);
}

const std::vector<double>& PredictionValues::get_values(std::size_t node) const {
  return values.at(node);
}

double PredictionValues::empty(std::size_t node) const {
  return values.at(node).empty();
}

const std::vector<std::vector<double>>& PredictionValues::get_all_values() const {
  return values;
}

const size_t PredictionValues::get_num_nodes() const {
  return num_nodes;
}

const size_t PredictionValues::get_num_types() const {
  return num_types;
}

void PredictionValues::clear() {
  num_nodes = 0;
  num_types = 0;
  values.clear();
}
