/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_PREDICTIONVALUES_H
#define GRF_PREDICTIONVALUES_H


#include <vector>
#include <string>
#include <map>

class PredictionValues {
public:
  PredictionValues();

  PredictionValues(const std::vector<std::vector<double>>& values,
                   size_t num_nodes,
                   size_t num_types);

  double empty(size_t node) const;

  double get(size_t node, size_t type) const;

  const std::vector<double>& get_values(size_t nodeID) const;

  const size_t get_num_nodes() const {
    return num_nodes;
  }

  const size_t get_num_types() const {
    return num_types;
  }

private:
  std::vector<std::vector<double>> values;
  size_t num_nodes;
  size_t num_types;
};


#endif //GRF_PREDICTIONVALUES_H
