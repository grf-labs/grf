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

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "utility.h"
#include "DefaultData.h"
#include "SparseData.h"

namespace grf {

void split_sequence(std::vector<uint>& result, uint start, uint end, uint num_parts) {

  result.reserve(num_parts + 1);

  // Return range if only 1 part
  if (num_parts == 1) {
    result.push_back(start);
    result.push_back(end + 1);
    return;
  }

  // Return vector from start to end+1 if more parts than elements
  if (num_parts > end - start + 1) {
    for (uint i = start; i <= end + 1; ++i) {
      result.push_back(i);
    }
    return;
  }

  uint length = (end - start + 1);
  uint part_length_short = length / num_parts;
  uint part_length_long = (uint) std::ceil(length / ((double) num_parts));
  uint cut_pos = length % num_parts;

  // Add long ranges
  for (uint i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
    result.push_back(i);
  }

  // Add short ranges
  for (uint i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
    result.push_back(i);
  }
}

bool equal_doubles(double first, double second, double epsilon) {
  if (std::isnan(first)) {
    return std::isnan(second);
  }
  return std::abs(first - second) < epsilon;
}

std::unique_ptr<Data> load_data(const std::string& file_name) {
  std::unique_ptr<Data> data = std::unique_ptr<Data>(new DefaultData());
  bool rounding_error = data->load_from_file(file_name);
  if (rounding_error) {
    throw std::runtime_error("A rounding error occurred while loading data from file.");
  }
  data->sort();
  return data;
}

std::unique_ptr<Data> load_sparse_data(const std::string& file_name) {
    std::unique_ptr<Data> data = std::unique_ptr<Data>(new SparseData());
  bool rounding_error = data->load_from_file(file_name);
  if (rounding_error) {
    throw std::runtime_error("A rounding error occurred while loading data from file.");
  }
  data->sort();
  return data;
}

} // namespace grf
