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

#include "DefaultData.h"
#include "utility.h"

DefaultData::DefaultData() :
    DefaultData(NULL, 0, 0) {}

DefaultData::DefaultData(double* data,
                         size_t num_rows,
                         size_t num_cols) :
    data(data) {
  this->num_rows = num_rows;
  this->num_cols = num_cols;
}

DefaultData::~DefaultData() {
  if (!external_data) {
    delete[] data;
  }
}

double DefaultData::get(size_t row, size_t col) const {
  return data[col * num_rows + row];
}

void DefaultData::reserve_memory() {
  data = new double[num_cols * num_rows];
}

void DefaultData::set(size_t col, size_t row, double value, bool& error) {
  data[col * num_rows + row] = value;
}
