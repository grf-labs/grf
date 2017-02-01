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

#ifndef GRADIENTFOREST_OBSERVATIONS_H
#define GRADIENTFOREST_OBSERVATIONS_H


#include <vector>
#include <string>
#include <map>

class Observations {
public:
  Observations(const std::map<std::string, std::vector<double>>& observations_by_type,
               size_t num_samples);

  const std::vector<double>& get(std::string type) const;

  const std::map<std::string, std::vector<double>>& get_observations_by_type() const {
    return observations_by_type;
  }

  size_t get_num_samples() const {
    return num_samples;
  }

  static const std::string OUTCOME;
  static const std::string TREATMENT;
  static const std::string INSTRUMENT;

private:
  std::map<std::string, std::vector<double>> observations_by_type;
  size_t num_samples;
};


#endif //GRADIENTFOREST_OBSERVATIONS_H
