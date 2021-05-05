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
#include <fstream>
#include <stdexcept>

#include "utility.h"

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

std::pair<std::vector<double>, std::vector<size_t>> load_data(const std::string& file_name) {
  size_t num_rows = 0;
  size_t num_cols = 0;

  // Open input file
  std::ifstream input_file;
  input_file.open(file_name);
  if (!input_file.good()) {
    throw std::runtime_error("Could not open input file.");
  }

  // Count number of rows
  size_t line_count = 0;
  std::string line;
  std::string first_line;
  while (getline(input_file, line)) {
    if (line_count == 0) {
      first_line = line;
    }
    ++line_count;
  }

  num_rows = line_count;
  input_file.close();
  input_file.open(file_name);

  // Read the first line to determine the number of columns.
  std::string dummy_token;
  std::stringstream first_line_stream(first_line);
  while (first_line_stream >> dummy_token) {
    num_cols++;
  }

  // Read the entire contents.
  std::vector<double> storage(num_rows * num_cols);
  line.clear();
  size_t row = 0;
  while (getline(input_file, line)) {
    std::string token;
    std::stringstream line_stream(line);
    size_t column = 0;
    while (line_stream >> token) {
      storage.at(column * num_rows + row) = std::stod(token);
      ++column;
    }
    if (column > num_cols) {
      throw std::runtime_error("Could not open input file. Too many columns in a row.");
    } else if (column < num_cols) {
      throw std::runtime_error("Could not open input file. Too few columns in a row. Are all values numeric?");
    }
    ++row;
  }
  num_rows = row;
  input_file.close();

  std::vector<size_t> dim {num_rows, num_cols};

  return std::make_pair(storage, dim);
}

void set_data(std::pair<std::vector<double>, std::vector<size_t>>& data, size_t row, size_t col, double value) {
  const std::vector<size_t>& dim = data.second;
  size_t num_rows = dim[0];

  data.first.at(col * num_rows + row) = value;
}

} // namespace grf
