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

#ifndef GRF_OBSERVATIONS_H
#define GRF_OBSERVATIONS_H


#include <vector>
#include <string>

class Observations {
public:
  Observations();

  Observations(const std::vector<std::vector<double>>& observations_by_type,
               size_t num_samples);

  double get(std::size_t type, size_t sample) const;
  const std::vector<std::vector<double>>& get_observations_by_type() const;

  size_t get_num_samples() const;
  size_t get_num_types() const;

  static const std::size_t OUTCOME;
  static const std::size_t TREATMENT;
  static const std::size_t INSTRUMENT;

private:
  std::vector<std::vector<double>> observations_by_type;
  size_t num_samples;
};


#endif //GRF_OBSERVATIONS_H
