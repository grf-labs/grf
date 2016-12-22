//
// Created by Julie Noelle Tibshirani on 12/22/16.
//

#ifndef GRADIENTFOREST_FILEUTILITIES_H
#define GRADIENTFOREST_FILEUTILITIES_H


#include <vector>
#include <string>

class FileTestUtilities {
public:
  static std::vector<std::vector<double>> readCsvFile(std::string file_name);


};


#endif //GRADIENTFOREST_FILEUTILITIES_H
