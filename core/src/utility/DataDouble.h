#ifndef GRADIENTFOREST_DATADOUBLE_H_
#define GRADIENTFOREST_DATADOUBLE_H_

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

#endif /* GRADIENTFOREST_DATADOUBLE_H_ */
