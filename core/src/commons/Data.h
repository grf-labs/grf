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

  size_t get_variable_id(std::string variable_name);

  void reserve_memory();
  void set(size_t col, size_t row, double value, bool& error);

  void add_sparse_data(unsigned char *sparse_data, size_t num_cols_sparse);

  bool load_from_file(std::string filename);
  bool load_from_whitespace_file(std::ifstream& input_file, std::string header_line);
  bool load_from_other_file(std::ifstream& input_file, std::string header_line, char seperator);

  void get_all_values(std::vector<double>& all_values, const std::vector<size_t>& samples, size_t var);

  size_t get_index(size_t row, size_t col) const {
    if (col < num_cols_no_sparse) {
      return index_data[col * num_rows + row];
    } else {
      // Get data out of sparse storage. -1 because of GenABEL coding.
      size_t idx = (col - num_cols_no_sparse) * num_rows_rounded + row;
      size_t result = (((sparse_data[idx / 4]&  mask[idx % 4]) >> offset[idx % 4]) - 1);

      // TODO: Better way to treat missing values?
      if (result > 2) {
        return 0;
      } else {
        return result;
      }
    }
  }

  double get_unique_data_value(size_t var, size_t index) const {
    if (var < num_cols_no_sparse) {
      return unique_data_values[var][index];
    } else {
      // For GWAS data the index is the value
      return (index);
    }
  }

  size_t get_num_unique_data_values(size_t var) const {
    if (var < num_cols_no_sparse) {
      return unique_data_values[var].size();
    } else {
      // For GWAS data 0,1,2
      return (3);
    }
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
    if (sparse_data == 0 || max_num_unique_values > 3) {
      // If no sparse data or one variable with more than 3 unique values, return that value
      return max_num_unique_values;
    } else {
      // If sparse data and no variable with more than 3 unique values, return 3
      return 3;
    }
  }

protected:
  std::vector<std::string> variable_names;
  size_t num_rows;
  size_t num_rows_rounded;
  size_t num_cols;

  unsigned char* sparse_data;
  size_t num_cols_no_sparse;

  bool externalData;

  size_t* index_data;
  std::vector<std::vector<double>> unique_data_values;
  size_t max_num_unique_values;

  double* data;

private:
  DISALLOW_COPY_AND_ASSIGN(Data);
};

#endif /* GRF_DATA_H_ */
