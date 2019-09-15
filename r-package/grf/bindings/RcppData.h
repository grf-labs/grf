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

#ifndef GRF_RCPPDATA_H
#define GRF_RCPPDATA_H

#include <Rcpp.h>

#include "commons/Data.h"
#include "commons/utility.h"

using namespace grf;

class RcppData final: public Data {
public:
  RcppData(Rcpp::NumericMatrix& data, size_t num_rows, size_t num_cols);

  double get(size_t row, size_t col) const;

  void reserve_memory();

  void set(size_t col, size_t row, double value, bool& error);

private:
  Rcpp::NumericMatrix data;

  DISALLOW_COPY_AND_ASSIGN(RcppData);
};

#endif //GRF_RCPPDATA_H
