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

#ifndef GRF_SPARSEDATA_H
#define GRF_SPARSEDATA_H


#include "Data.h"
#include "Eigen/Sparse"

class SparseData: public Data {
public:
  SparseData();

  SparseData(Eigen::SparseMatrix<double>* data,
             std::vector<std::string> variable_names,
             size_t num_rows,
             size_t num_cols);

  virtual ~SparseData();

  double get(size_t row, size_t col) const;

  void reserve_memory();
  void set(size_t col, size_t row, double value, bool& error);

protected:
  Eigen::SparseMatrix<double>* data;

private:
  DISALLOW_COPY_AND_ASSIGN(SparseData);
};


#endif //GRF_SPARSEDATA_H
