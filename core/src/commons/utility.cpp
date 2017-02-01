
#include <math.h>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "utility.h"
#include "globals.h"
#include "Data.h"

void equalSplit(std::vector<uint>& result, uint start, uint end, uint num_parts) {

  result.reserve(num_parts + 1);

  // Return range if only 1 part
  if (num_parts == 1) {
    result.push_back(start);
    result.push_back(end + 1);
    return;
  }

  // Return vector from start to end+1 if more parts than elements
  if (num_parts > end - start + 1) {
    for (uint i = start; i <= end + 1; ++i) {
      result.push_back(i);
    }
    return;
  }

  uint length = (end - start + 1);
  uint part_length_short = length / num_parts;
  uint part_length_long = (uint) ceil(length / ((double) num_parts));
  uint cut_pos = length % num_parts;

  // Add long ranges
  for (uint i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
    result.push_back(i);
  }

  // Add short ranges
  for (uint i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
    result.push_back(i);
  }
}

void loadDoubleVectorFromFile(std::vector<double>& result, std::string filename) {

  // Open input file
  std::ifstream input_file;
  input_file.open(filename);
  if (!input_file.good()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  // Read the first line, ignore the rest
  std::string line;
  getline(input_file, line);
  std::stringstream line_stream(line);
  double token;
  while (line_stream >> token) {
    result.push_back(token);
  }
  input_file.close();
}

std::string uintToString(uint number) {
  return std::to_string(number);
}

std::string beautifyTime(uint seconds) {
  std::string result;

  // Add seconds, minutes, hours, days if larger than zero
  uint out_seconds = (uint) seconds % 60;
  result = uintToString(out_seconds) + " seconds";
  uint out_minutes = (seconds / 60) % 60;
  if (seconds / 60 == 0) {
    return result;
  } else if (out_minutes == 1) {
    result = "1 minute, " + result;
  } else {
    result = uintToString(out_minutes) + " minutes, " + result;
  }
  uint out_hours = (seconds / 3600) % 24;
  if (seconds / 3600 == 0) {
    return result;
  } else if (out_hours == 1) {
    result = "1 hour, " + result;
  } else {
    result = uintToString(out_hours) + " hours, " + result;
  }
  uint out_days = (seconds / 86400);
  if (out_days == 0) {
    return result;
  } else if (out_days == 1) {
    result = "1 day, " + result;
  } else {
    result = uintToString(out_days) + " days, " + result;
  }
  return result;
}

size_t roundToNextMultiple(size_t value, uint multiple) {

  if (multiple == 0) {
    return value;
  }

  int remainder = value % multiple;
  if (remainder == 0) {
    return value;
  }

  return value + multiple - remainder;
}

void splitString(std::vector<std::string>& result, std::string input, char split_char) {

  std::istringstream ss(input);
  std::string token;

  while (std::getline(ss, token, split_char)) {
    result.push_back(token);
  }
}

bool equalDoubles(double first, double second, double epsilon) {
  if (std::isnan(first)) {
    return std::isnan(second);
  }
  return std::abs(first - second) < epsilon;
}

Data* loadDataFromFile(std::string file_name) {
  Data* data = new Data();

  bool rounding_error = data->loadFromFile(file_name);
  if (rounding_error) {
    std::cerr << "A rounding error occurred while loading data from file." << std::endl;
  }

  data->sort();
  return data;
}