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

#include "commons/Data.h"
#include "forest/Forest.h"

Forest Forest::create(std::vector<std::shared_ptr<Tree>> trees,
                      Data* data,
                      std::unordered_map<size_t, size_t> observables) {
  size_t num_types = observables.size();
  size_t num_samples = data->get_num_rows();

  std::vector<std::vector<double>> observations_by_type(num_types);
  for (auto it : observables) {
    size_t type = it.first;
    size_t index = it.second;

    observations_by_type[type].resize(num_samples);
    for (size_t row = 0; row < num_samples; ++row) {
      observations_by_type[type][row] = data->get(row, index);
    }
  }

  Observations observations(observations_by_type, num_samples);
  return Forest(trees, observations, data->get_num_cols());
}

Forest::Forest(const std::vector<std::shared_ptr<Tree>>& trees,
               const Observations& observations,
               size_t num_variables):
  trees(trees),
  observations(observations),
  num_variables(num_variables) {}
