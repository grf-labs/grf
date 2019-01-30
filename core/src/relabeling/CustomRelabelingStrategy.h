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

#ifndef GRF_CUSTOMRELABELINGSTRATEGY_H
#define GRF_CUSTOMRELABELINGSTRATEGY_H

#include <numeric>
 
#include "RelabelingStrategy.h"

class CustomRelabelingStrategy: public RelabelingStrategy {
public:
  std::unordered_map<size_t, double> relabel(
      const std::vector<size_t>& samples,
      const Data* data);
};

std::vector<double> logrankScores(const std::vector<double>& time, const std::vector<double>& status);

template<typename T>
std::vector<size_t> order(const std::vector<T>& values, bool decreasing) {
  // Create index vector
  std::vector<size_t> indices(values.size());
  std::iota(indices.begin(), indices.end(), 0);
  
  // Sort index vector based on value vector
  if (decreasing) {
    std::sort(std::begin(indices), std::end(indices), [&](size_t i1, size_t i2) {return values[i1] > values[i2];});
  } else {
    std::sort(std::begin(indices), std::end(indices), [&](size_t i1, size_t i2) {return values[i1] < values[i2];});
  }
  return indices;
}


#endif //GRF_CUSTOMRELABELINGSTRATEGY_H
