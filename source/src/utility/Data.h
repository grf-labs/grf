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

protected:
  std::vector<std::string> variable_names;
  size_t num_rows;
  size_t num_rows_rounded;
  size_t num_cols;

  unsigned char* sparse_data;
  size_t num_cols_no_sparse;

  bool externalData;

  // TODO: Need normal data still?
  // TODO: Dont need all three!
  size_t* index_data;
  // TODO: Really need this? Just use index everywhere!
  std::vector<std::vector<double>> unique_data_values;

private:
  DISALLOW_COPY_AND_ASSIGN(Data);
};

#endif /* DATA_H_ */
