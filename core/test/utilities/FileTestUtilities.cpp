#include <fstream>
#include "FileTestUtilities.h"

std::vector<std::vector<double>> FileTestUtilities::readCsvFile(std::string file_name) {
  std::ifstream file = std::ifstream();
  file.open(file_name, std::ios::binary);

  std::string delimiter = " ";
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

void FileTestUtilities::writeCsvFile(std::string file_name, std::vector<std::vector<double>> contents) {
  std::ofstream file;
  file.open(file_name, std::ios::binary);

  for (auto &line : contents) {
    for (int i = 0; i <line.size(); i++) {
      if (i > 0) {
        file << ", ";
      }
      file << line[i];
    }
    file << std::endl;
  }
}

