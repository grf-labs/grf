/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include <algorithm>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <sstream>

#include "Data.h"

namespace grf {

Data::Data() :
    num_rows(0),
    num_cols(0),
    index_data(),
    max_num_unique_values(0),
    outcome_index(),
    treatment_index(),
    instrument_index(),
    weight_index() {}

bool Data::load_from_file(const std::string& filename) {
  bool result;

  // Open input file
  std::ifstream input_file;
  input_file.open(filename);
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
  input_file.open(filename);

  // Find out if comma, semicolon or whitespace seperated and call appropriate method
  if (first_line.find(',') != std::string::npos) {
    result = load_from_other_file(input_file, first_line, ',');
  } else if (first_line.find(';') != std::string::npos) {
    result = load_from_other_file(input_file, first_line, ';');
  } else {
    result = load_from_whitespace_file(input_file, first_line);
  }

  input_file.close();
  return result;
}

bool Data::load_from_whitespace_file(std::ifstream& input_file,
                                     const std::string& first_line) {
  // Read the first line to determine the number of columns.
  std::string dummy_token;
  std::stringstream first_line_stream(first_line);
  while (first_line_stream >> dummy_token) {
    num_cols++;
  }

  // Read the entire contents.
  reserve_memory();
  bool error = false;
  std::string line;
  size_t row = 0;
  while (getline(input_file, line)) {
    double token;
    std::stringstream line_stream(line);
    size_t column = 0;
    while (line_stream >> token) {
      set(column, row, token, error);
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
  return error;
}

bool Data::load_from_other_file(std::ifstream& input_file,
                                const std::string& first_line,
                                char seperator) {
  // Read the first line to determine the number of columns.
  std::string dummy_token;
  std::stringstream first_line_stream(first_line);
  while (getline(first_line_stream, dummy_token, seperator)) {
    num_cols++;
  }

  // Read the entire contents.
  reserve_memory();
  bool error = false;
  std::string line;
  size_t row = 0;
  while (getline(input_file, line)) {
    std::string token_string;
    double token;
    std::stringstream line_stream(line);
    size_t column = 0;
    while (getline(line_stream, token_string, seperator)) {
      std::stringstream token_stream(token_string);
      token_stream >> token;
      set(column, row, token, error);
      ++column;
    }
    ++row;
  }
  num_rows = row;
  return error;
}

void Data::set_outcome_index(size_t index) {
  this->outcome_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_treatment_index(size_t index) {
  this->treatment_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_instrument_index(size_t index) {
  this->instrument_index = index;
  disallowed_split_variables.insert(index);
}

void Data::set_weight_index(size_t index) {
  this->weight_index = index;
  disallowed_split_variables.insert(index);
}

void Data::get_all_values(std::vector<double>& all_values, const std::vector<size_t>& samples, size_t var) const {
  all_values.resize(samples.size());
  for (size_t i = 0; i < samples.size(); i++) {
    size_t sample = samples[i];
    all_values[i] = get(sample, var);
  }
  std::sort(all_values.begin(), all_values.end());
  all_values.erase(unique(all_values.begin(), all_values.end()), all_values.end());
}

size_t Data::get_index(size_t row, size_t col) const {
  return index_data[col * num_rows + row];
}

void Data::sort() {
  // Reserve memory
  index_data.resize(num_cols * num_rows);

  // For all columns, get unique values and save index for each observation
  for (size_t col = 0; col < num_cols; ++col) {

    // Get all unique values
    std::vector<double> unique_values(num_rows);
    for (size_t row = 0; row < num_rows; ++row) {
      unique_values[row] = get(row, col);
    }
    std::sort(unique_values.begin(), unique_values.end());
    unique_values.erase(unique(unique_values.begin(), unique_values.end()), unique_values.end());

    // Get index of unique value
    for (size_t row = 0; row < num_rows; ++row) {
      size_t idx =
          std::lower_bound(unique_values.begin(), unique_values.end(), get(row, col)) - unique_values.begin();
      index_data[col * num_rows + row] = idx;
    }

    // Save unique values
    unique_data_values.push_back(unique_values);
    if (unique_values.size() > max_num_unique_values) {
      max_num_unique_values = unique_values.size();
    }
  }
}

double Data::get_unique_data_value(size_t var, size_t index) const {
  return unique_data_values[var][index];
}

size_t Data::get_num_unique_data_values(size_t var) const {
  return unique_data_values[var].size();
}

size_t Data::get_num_cols() const {
  return num_cols;
}

size_t Data::get_num_rows() const {
  return num_rows;
}

size_t Data::get_max_num_unique_values() const {
  return max_num_unique_values;
}

double Data::get_outcome(size_t row) const {
  return get(row, outcome_index.value());
}

double Data::get_treatment(size_t row) const {
  return get(row, treatment_index.value());
}

double Data::get_instrument(size_t row) const {
  return get(row, instrument_index.value());
}

double Data::get_weight(size_t row) const {
  if (weight_index.has_value()) {
    return get(row, weight_index.value());
  } else {
    return 1.0;
  }
}

const std::set<size_t>& Data::get_disallowed_split_variables() const {
  return disallowed_split_variables;
}

} // namespace grf
