/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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
#include "FileTestUtilities.h"

std::vector<std::vector<double>> FileTestUtilities::read_csv_file(const std::string& file_name) {
  std::ifstream file;
  file.open(file_name, std::ios::binary);

  std::string delimiter = ",";
  std::vector<std::vector<double>> result;
  std::string line;

  while (std::getline(file, line)) {
    size_t position = 0;
    std::string token;

    std::vector<double> split_line;
    while ((position = line.find(delimiter)) != std::string::npos) {
      token = line.substr(0, position);
      split_line.push_back(std::stod(token));
      line.erase(0, position + delimiter.length());
    }

    if (position != line.length()) {
      token = line.substr(0, line.length());
      split_line.push_back(std::stod(token));
    }

  result.push_back(split_line);
  }

  file.close();
  return result;
}

void FileTestUtilities::write_csv_file(const std::string& file_name,
                                       const std::vector<std::vector<double>>& contents) {
  std::ofstream file;
  file.open(file_name, std::ios::binary);

  for (auto& line : contents) {
    for (size_t i = 0; i < line.size(); ++i) {
      if (i > 0) {
        file << ", ";
      }
      file << line[i];
    }
    file << std::endl;
  }
  file.close();
}
