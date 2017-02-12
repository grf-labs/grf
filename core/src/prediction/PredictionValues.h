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

#ifndef GRADIENTFOREST_PREDICTIONVALUES_H
#define GRADIENTFOREST_PREDICTIONVALUES_H


#include <vector>
#include <string>
#include <map>

class PredictionValues {
public:
  PredictionValues();

  PredictionValues(const std::map<std::string, std::vector<double>>& values_by_type,
                   size_t num_nodes);

  const std::vector<double>& get(std::string type) const;

  const size_t get_num_nodes() const {
    return num_nodes;
  }

  const std::map<std::string, std::vector<double>>& get_values_by_type() const {
    return values_by_type;
  }

  const bool empty() const {
    return num_nodes == 0;
  }

private:
  std::map<std::string, std::vector<double>> values_by_type;
  size_t num_nodes;
};


#endif //GRADIENTFOREST_PREDICTIONVALUES_H
