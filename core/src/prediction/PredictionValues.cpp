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

#include "PredictionValues.h"

PredictionValues::PredictionValues():
  num_nodes(0),
  num_types(0) {}

PredictionValues::PredictionValues(const std::vector<std::vector<double>>& values,
                                   size_t num_nodes,
                                   size_t num_types):
  values(values),
  num_nodes(num_nodes),
  num_types(num_types) {}

double PredictionValues::empty(std::size_t node) const {
  return values.at(node).empty();
}

double PredictionValues::get(std::size_t node, size_t type) const {
  return values.at(node).at(type);
}


