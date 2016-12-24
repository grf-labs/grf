#ifndef GRADIENTFOREST_DATA_H_
#define GRADIENTFOREST_DATA_H_

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

  size_t getVariableID(std::string variable_name);

  void reserveMemory();
  void set(size_t col, size_t row, double value, bool& error);

  void addSparseData(unsigned char* sparse_data, size_t num_cols_sparse);

  bool loadFromFile(std::string filename);
  bool loadFromFileWhitespace(std::ifstream& input_file, std::string header_line);
  bool loadFromFileOther(std::ifstream& input_file, std::string header_line, char seperator);

  void getAllValues(std::vector<double>& all_values, std::vector<size_t>& sampleIDs, size_t varID);

  size_t getIndex(size_t row, size_t col) const {
    if (col < num_cols_no_sparse) {
      return index_data[col * num_rows + row];
    } else {
      // Get data out of sparse storage. -1 because of GenABEL coding.
      size_t idx = (col - num_cols_no_sparse) * num_rows_rounded + row;
      size_t result = (((sparse_data[idx / 4] & mask[idx % 4]) >> offset[idx % 4]) - 1);

      // TODO: Better way to treat missing values?
      if (result > 2) {
        return 0;
      } else {
        return result;
      }
    }
  }

  double getUniqueDataValue(size_t varID, size_t index) const {
    if (varID < num_cols_no_sparse) {
      return unique_data_values[varID][index];
    } else {
      // For GWAS data the index is the value
      return (index);
    }
  }

  size_t getNumUniqueDataValues(size_t varID) const {
    if (varID < num_cols_no_sparse) {
      return unique_data_values[varID].size();
    } else {
      // For GWAS data 0,1,2
      return (3);
    }
  }

  void sort();

  const std::vector<std::string>& getVariableNames() const {
    return variable_names;
  }
  size_t getNumCols() const {
    return num_cols;
  }
  size_t getNumRows() const {
    return num_rows;
  }

  size_t getMaxNumUniqueValues() const {
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

#endif /* GRADIENTFOREST_DATA_H_ */
