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

#ifndef GRF_DATA_H_
#define GRF_DATA_H_

#include <vector>
#include <iostream>

#include "globals.h"

class Data {
public:
  Data();

  Data(double* data,
       std::vector<std::string> variable_names,
       size_t num_rows,
       size_t num_cols);

  ~Data();

  double get(size_t row, size_t col) const;

  void reserve_memory();
  void set(size_t col, size_t row, double value, bool& error);

  bool load_from_file(std::string filename);
  bool load_from_whitespace_file(std::ifstream& input_file, std::string header_line);
  bool load_from_other_file(std::ifstream& input_file, std::string header_line, char seperator);

  void get_all_values(std::vector<double>& all_values, const std::vector<size_t>& samples, size_t var);

  size_t get_index(size_t row, size_t col) const {
    return index_data[col * num_rows + row];
  }

  double get_unique_data_value(size_t var, size_t index) const {
    return unique_data_values[var][index];
  }

  size_t get_num_unique_data_values(size_t var) const {
    return unique_data_values[var].size();
  }

  void sort();

  const std::vector<std::string>& get_variable_names() const {
    return variable_names;
  }
  size_t get_num_cols() const {
    return num_cols;
  }
  size_t get_num_rows() const {
    return num_rows;
  }

  size_t get_max_num_unique_values() const {
    return max_num_unique_values;
  }

protected:
  std::vector<std::string> variable_names;
  size_t num_rows;
  size_t num_rows_rounded;
  size_t num_cols;

  bool externalData;

  size_t* index_data;
  std::vector<std::vector<double>> unique_data_values;
  size_t max_num_unique_values;

  double* data;

private:
  DISALLOW_COPY_AND_ASSIGN(Data);
};

#endif /* GRF_DATA_H_ */
