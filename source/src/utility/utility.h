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

#ifndef UTILITY_H_
#define UTILITY_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <unordered_map>

#include "globals.h"
#include "Data.h"

/**
 * Split sequence start..end in num_parts parts with sizes as equal as possible.
 * @param result Result vector of size num_parts+1. Ranges for the parts are then result[0]..result[1]-1, result[1]..result[2]-1, ..
 * @param start minimum value
 * @param end maximum value
 * @param num_parts number of parts
 */
void equalSplit(std::vector<uint>& result, uint start, uint end, uint num_parts);

/**
 * Write a 1d vector to filestream. First the size is written as size_t, then all vector elements.
 * @param vector Vector with elements of type T to write to file.
 * @param file ofstream object to write to.
 */

/**
 * Write a 1d vector to filestream. First the size is written, then all vector elements.
 * @param vector Vector of type T to save
 * @param file ofstream object to write to.
 */
template<typename T>
inline void saveVector1D(std::vector<T>& vector, std::ofstream& file) {
  // Save length
  size_t length = vector.size();
  file.write((char*) &length, sizeof(length));
  file.write((char*) vector.data(), length * sizeof(T));
}

template<>
inline void saveVector1D(std::vector<bool>& vector, std::ofstream& file) {
  // Save length
  size_t length = vector.size();
  file.write((char*) &length, sizeof(length));

  // Save vector
  for (size_t i = 0; i < vector.size(); ++i) {
    bool v = vector[i];
    file.write((char*) &v, sizeof(v));
  }
}

/**
 * Read a 1d vector written by saveVector1D() from filestream.
 * @param result Result vector with elements of type T.
 * @param file ifstream object to read from.
 */
template<typename T>
inline void readVector1D(std::vector<T>& result, std::ifstream& file) {
  // Read length
  size_t length;
  file.read((char*) &length, sizeof(length));
  result.resize(length);
  file.read((char*) result.data(), length * sizeof(T));
}

template<>
inline void readVector1D(std::vector<bool>& result, std::ifstream& file) {
  // Read length
  size_t length;
  file.read((char*) &length, sizeof(length));

  // Read vector.
  for (size_t i = 0; i < length; ++i) {
    bool temp;
    file.read((char*) &temp, sizeof(temp));
    result.push_back(temp);
  }
}

/**
 * Write a 2d vector to filestream. First the size of the first dim is written as size_t, then for all inner vectors the size and elements.
 * @param vector Vector of vectors of type T to write to file.
 * @param file ofstream object to write to.
 */
template<typename T>
inline void saveVector2D(std::vector<std::vector<T>>& vector, std::ofstream& file) {
  // Save length of first dim
  size_t length = vector.size();
  file.write((char*) &length, sizeof(length));

  // Save outer vector
  for (auto& inner_vector : vector) {
    // Save inner vector
    saveVector1D(inner_vector, file);
  }
}

/**
 * Read a 2d vector written by saveVector2D() from filestream.
 * @param result Result vector of vectors with elements of type T.
 * @param file ifstream object to read from.
 */
template<typename T>
inline void readVector2D(std::vector<std::vector<T>>& result, std::ifstream& file) {
  // Read length of first dim
  size_t length;
  file.read((char*) &length, sizeof(length));
  result.resize(length);

  // Read outer vector
  for (size_t i = 0; i < length; ++i) {
    // Read inner vector
    readVector1D(result[i], file);
  }
}

/**
 * Read a double vector from text file. Reads only the first line.
 * @param result Result vector of doubles with contents
 * @param filename filename of input file
 */
void loadDoubleVectorFromFile(std::vector<double>& result, std::string filename);

/**
 * Convert a unsigned integer to string
 * @param number Number to convert
 * @return Converted number as string
 */
std::string uintToString(uint number);

/**
 * Beautify output of time.
 * @param seconds Time in seconds
 * @return Time in days, hours, minutes and seconds as string
 */
std::string beautifyTime(uint seconds);

/**
 * Round up to next multiple of a number.
 * @param value Value to be rounded up.
 * @param multiple Number to multiply.
 * @return Rounded number
 */
size_t roundToNextMultiple(size_t value, uint multiple);

/**
 * Split string in parts separated by character.
 * @param result Splitted string
 * @param input String to be splitted
 * @param split_char Char to separate parts
 */
void splitString(std::vector<std::string>& result, std::string input, char split_char);


/**
 * Check if all values in double vector are positive integers.
 * @param all_values Double vector to check
 * @return True if all values are positive integers
 */
bool checkPositiveIntegers(std::vector<double>& all_values);

#endif /* UTILITY_H_ */
