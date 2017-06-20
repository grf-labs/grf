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

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>

#include "Data.h"
#include "utility.h"

Data::Data() :
  Data(NULL, std::vector<std::string>(), 0, 0) {}

Data::Data(double* data,
           std::vector<std::string> variable_names,
           size_t num_rows,
           size_t num_cols):
    variable_names(variable_names),
    num_rows(num_rows),
    num_rows_rounded(0),
    num_cols(num_cols),
    sparse_data(0),
    num_cols_no_sparse(num_cols),
    externalData(true),
    index_data(0),
    max_num_unique_values(0),
    data(data) {}

Data::~Data() {
  if (index_data != 0) {
    delete[] index_data;
  }
  if (!externalData) {
    delete[] data;
  }
}

size_t Data::get_variable_id(std::string variable_name) {
  std::vector<std::string>::iterator it = std::find(variable_names.begin(), variable_names.end(), variable_name);
  if (it == variable_names.end()) {
    throw std::runtime_error("Variable " + variable_name + " not found.");
  }
  return (std::distance(variable_names.begin(), it));
}

void Data::add_sparse_data(unsigned char *sparse_data, size_t num_cols_sparse) {
  num_cols = num_cols_no_sparse + num_cols_sparse;
  num_rows_rounded = round_to_next_multiple(num_rows, 4);
  this->sparse_data = sparse_data;
}

bool Data::load_from_file(std::string filename) {

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
  while (getline(input_file, line)) {
    ++line_count;
  }
  num_rows = line_count - 1;
  input_file.close();
  input_file.open(filename);

  // Check if comma, semicolon or whitespace seperated
  std::string header_line;
  getline(input_file, header_line);

  // Find out if comma, semicolon or whitespace seperated and call appropriate method
  if (header_line.find(",") != std::string::npos) {
    result = load_from_other_file(input_file, header_line, ',');
  } else if (header_line.find(";") != std::string::npos) {
    result = load_from_other_file(input_file, header_line, ';');
  } else {
    result = load_from_whitespace_file(input_file, header_line);
  }

  externalData = false;
  input_file.close();
  return result;
}

bool Data::load_from_whitespace_file(std::ifstream& input_file, std::string header_line) {

  // Read header
  std::string header_token;
  std::stringstream header_line_stream(header_line);
  while (header_line_stream >> header_token) {
    variable_names.push_back(header_token);
  }
  num_cols = variable_names.size();
  num_cols_no_sparse = num_cols;

  // Read body
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

bool Data::load_from_other_file(std::ifstream& input_file, std::string header_line, char seperator) {

  // Read header
  std::string header_token;
  std::stringstream header_line_stream(header_line);
  while (getline(header_line_stream, header_token, seperator)) {
    variable_names.push_back(header_token);
  }
  num_cols = variable_names.size();
  num_cols_no_sparse = num_cols;

  // Read body
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

void Data::get_all_values(std::vector<double>& all_values, const std::vector<size_t>& samples, size_t var) {

  // All values for var (no duplicates) for given samples
  if (var < num_cols_no_sparse) {

    all_values.reserve(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
      all_values.push_back(get(samples[i], var));
    }
    std::sort(all_values.begin(), all_values.end());
    all_values.erase(unique(all_values.begin(), all_values.end()), all_values.end());
  } else {
    // If GWA data just use 0, 1, 2
    all_values = std::vector<double>( { 0, 1, 2 });
  }

}

double Data::get(size_t row, size_t col) const {
  if (col < num_cols_no_sparse) {
    return data[col * num_rows + row];
  } else {
    // Get data out of sparse storage. -1 because of GenABEL coding.
    size_t idx = (col - num_cols_no_sparse) * num_rows_rounded + row;
    double result = (((sparse_data[idx / 4]&  mask[idx % 4]) >> offset[idx % 4]) - 1);
    return result;
  }
}

void Data::sort() {

  // Reserve memory
  index_data = new size_t[num_cols_no_sparse * num_rows];

  // For all columns, get unique values and save index for each observation
  for (size_t col = 0; col < num_cols_no_sparse; ++col) {

    // Get all unique values
    std::vector<double> unique_values(num_rows);
    for (size_t row = 0; row < num_rows; ++row) {
      unique_values[row] = get(row, col);
    }
    std::sort(unique_values.begin(), unique_values.end());
    unique_values.erase(unique(unique_values.begin(), unique_values.end()), unique_values.end());

    // Get index of unique value
    for (size_t row = 0; row < num_rows; ++row) {
      size_t idx = std::lower_bound(unique_values.begin(), unique_values.end(), get(row, col)) - unique_values.begin();
      index_data[col * num_rows + row] = idx;
    }

    // Save unique values
    unique_data_values.push_back(unique_values);
    if (unique_values.size() > max_num_unique_values) {
      max_num_unique_values = unique_values.size();
    }
  }
}

void Data::reserve_memory() {
  data = new double[num_cols * num_rows];
}

void Data::set(size_t col, size_t row, double value, bool& error) {
  data[col * num_rows + row] = value;
}
