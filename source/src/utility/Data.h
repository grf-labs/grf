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

#ifndef DATA_H_
#define DATA_H_

#include <vector>
#include <iostream>

#include "globals.h"

class Data {
public:
  Data();
  virtual ~Data();

  virtual double get(size_t row, size_t col) const = 0;

  size_t getVariableID(std::string variable_name);

  virtual void reserveMemory() = 0;
  virtual void set(size_t col, size_t row, double value, bool& error) = 0;

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

private:
  DISALLOW_COPY_AND_ASSIGN(Data);
};

#endif /* DATA_H_ */
