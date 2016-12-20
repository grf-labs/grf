/*-------------------------------------------------------------------------------
This file is part of Ranger.
    
Ranger is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Ranger is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Ranger. If not, see <http://www.gnu.org/licenses/>.

Written by: 

Marvin N. Wright
Institut für Medizinische Biometrie und Statistik
Universität zu Lübeck
Ratzeburger Allee 160
23562 Lübeck 

http://www.imbs-luebeck.de
wright@imbs.uni-luebeck.de
#-------------------------------------------------------------------------------*/

#ifndef DATADOUBLE_H_
#define DATADOUBLE_H_

#include "globals.h"
#include "utility.h"
#include "Data.h"

class DataDouble: public Data {
public:
  DataDouble();
  DataDouble(double* data, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols) :
      data(data) {
    this->variable_names = variable_names;
    this->num_rows = num_rows;
    this->num_cols = num_cols;
    this->num_cols_no_sparse = num_cols;
  }
  virtual ~DataDouble();

  double get(size_t row, size_t col) const {
    if (col < num_cols_no_sparse) {
      return data[col * num_rows + row];
    } else {
      // Get data out of sparse storage. -1 because of GenABEL coding.
      size_t idx = (col - num_cols_no_sparse) * num_rows_rounded + row;
      double result = (((sparse_data[idx / 4] & mask[idx % 4]) >> offset[idx % 4]) - 1);
      return result;
    }
  }

  void reserveMemory() {
    data = new double[num_cols * num_rows];
  }

  void set(size_t col, size_t row, double value, bool& error) {
    data[col * num_rows + row] = value;
  }

private:
  double* data;

  DISALLOW_COPY_AND_ASSIGN(DataDouble);
};

#endif /* DATADOUBLE_H_ */
