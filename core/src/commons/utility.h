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

#ifndef GRF_UTILITY_H_
#define GRF_UTILITY_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "Data.h"
#include "globals.h"

namespace grf {

/**
 * Split sequence start..end in num_parts parts with sizes as equal as possible.
 * @param result Result vector of size num_parts+1. Ranges for the parts are then result[0]..result[1]-1, result[1]..result[2]-1, ..
 * @param start minimum value
 * @param end maximum value
 * @param num_parts number of parts
 */
void split_sequence(std::vector<uint>& result, uint start, uint end, uint num_parts);

bool equal_doubles(double first, double second, double epsilon);

std::unique_ptr<Data> load_data(const std::string& file_name);

std::unique_ptr<Data> load_sparse_data(const std::string& file_name);

} // namespace grf

#endif /* GRF_UTILITY_H_ */
