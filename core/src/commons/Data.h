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
  virtual ~Data() {};

  virtual double get(size_t row, size_t col) const = 0;

  virtual void reserve_memory() = 0;
  virtual void set(size_t col, size_t row, double value, bool& error) = 0;

  virtual bool load_from_file(std::string filename) = 0;
  virtual bool load_from_whitespace_file(std::ifstream& input_file, std::string header_line) = 0;
  virtual bool load_from_other_file(std::ifstream& input_file, std::string header_line, char seperator)  = 0;

  virtual void get_all_values(std::vector<double>& all_values, const std::vector<size_t>& samples, size_t var)  = 0;

  virtual size_t get_index(size_t row, size_t col) const = 0;

  virtual double get_unique_data_value(size_t var, size_t index) const = 0;

  virtual size_t get_num_unique_data_values(size_t var) const = 0;

  virtual void sort() = 0;

  virtual const std::vector<std::string>& get_variable_names() const  = 0;

  virtual size_t get_num_cols() const = 0;
  virtual size_t get_num_rows() const = 0;

  virtual size_t get_max_num_unique_values() const = 0;
};

#endif /* GRF_DATA_H_ */
