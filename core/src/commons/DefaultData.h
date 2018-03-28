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

#ifndef GRF_DEFAULTDATA_H_
#define GRF_DEFAULTDATA_H_

#include <vector>
#include <iostream>

#include "Data.h"
#include "globals.h"

class DefaultData: public Data {
public:
  DefaultData();

  DefaultData(double* data,
       size_t num_rows,
       size_t num_cols);

  virtual ~DefaultData();

  double get(size_t row, size_t col) const;

  void reserve_memory();
  void set(size_t col, size_t row, double value, bool& error);

protected:
  double* data;

private:
  DISALLOW_COPY_AND_ASSIGN(DefaultData);
};

#endif /* GRF_DEFAULTDATA_H_ */
