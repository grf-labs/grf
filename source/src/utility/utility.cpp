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
}

// TODO: Remove
//void drawWithoutReplacementSkip(std::unordered_set<size_t>& result, std::mt19937_64& random_number_generator,
//    size_t max, std::vector<size_t>& skip, size_t num_samples) {
//
//  std::uniform_int_distribution<size_t> unif_dist(0, max - 1 - skip.size());
//  for (size_t i = 0; i < num_samples; ++i) {
//    size_t draw;
//    do {
//      draw = unif_dist(random_number_generator);
//      for (auto& skip_value : skip) {
//        if (draw >= skip_value) {
//          ++draw;
//        }
//      }
//    } while (result.count(draw) > 0);
//    result.insert(draw);
//  }
//}

void drawWithoutReplacementSkip(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    std::vector<size_t>& skip, size_t num_samples) {
  if (num_samples < max / 2) {
    drawWithoutReplacementSimple(result, random_number_generator, max, skip, num_samples);
  } else {
    drawWithoutReplacementKnuth(result, random_number_generator, max, skip, num_samples);
  }
}

void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    std::vector<size_t>& skip, size_t num_samples) {

  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(max, false);

  std::uniform_int_distribution<size_t> unif_dist(0, max - 1 - skip.size());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = unif_dist(random_number_generator);
      for (auto& skip_value : skip) {
        if (draw >= skip_value) {
          ++draw;
        }
      }
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

void drawWithoutReplacementKnuth(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    std::vector<size_t>& skip, size_t num_samples) {

  size_t size_no_skip = max - skip.size();
  result.resize(num_samples);
  double u;
  size_t final_value;

  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  size_t i = 0;
  size_t j = 0;
  while (i < num_samples) {
    u = distribution(random_number_generator);

    if ((size_no_skip - j) * u >= num_samples - i) {
      j++;
    } else {
      final_value = j;
      for (auto& skip_value : skip) {
        if (final_value >= skip_value) {
          ++final_value;
        }
      }
      result[i] = final_value;
      j++;
      i++;
    }
  }
}

void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    std::vector<size_t>& indizes, size_t num_samples, std::vector<double>& weights) {

  std::unordered_set<size_t> temp;
  std::discrete_distribution<> weighted_dist(weights.begin(), weights.end());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = weighted_dist(random_number_generator);
    } while (temp.count(indizes[draw]) > 0);
    temp.insert(indizes[draw]);
  }
  result.resize(num_samples);
  std::copy(temp.begin(), temp.end(), result.begin());
}

double mostFrequentValue(std::unordered_map<double, size_t>& class_count, std::mt19937_64 random_number_generator) {
  std::vector<double> major_classes;

  // Find maximum count
  size_t max_count = 0;
  for (auto& class_value : class_count) {
    if (class_value.second > max_count) {
      max_count = class_value.second;
      major_classes.clear();
      major_classes.push_back(class_value.first);
    } else if (class_value.second == max_count) {
      major_classes.push_back(class_value.first);
    }
  }

  if (major_classes.size() == 1) {
    return major_classes[0];
  } else {
    // Choose randomly
    std::uniform_int_distribution<size_t> unif_dist(0, major_classes.size() - 1);
    return major_classes[unif_dist(random_number_generator)];
  }
}

double computeConcordanceIndex(Data* data, std::vector<double>& sum_chf, size_t dependent_varID, size_t status_varID,
    std::vector<size_t>& sample_IDs) {

  // Compute concordance index
  double concordance = 0;
  double permissible = 0;
  for (size_t i = 0; i < sum_chf.size(); ++i) {
    size_t sample_i = i;
    if (!sample_IDs.empty()) {
      sample_i = sample_IDs[i];
    }
    double time_i = data->get(sample_i, dependent_varID);
    double status_i = data->get(sample_i, status_varID);

    for (size_t j = i + 1; j < sum_chf.size(); ++j) {
      size_t sample_j = j;
      if (!sample_IDs.empty()) {
        sample_j = sample_IDs[j];
      }
      double time_j = data->get(sample_j, dependent_varID);
      double status_j = data->get(sample_j, status_varID);

      if (time_i < time_j && status_i == 0) {
        continue;
      }
      if (time_j < time_i && status_j == 0) {
        continue;
      }
      if (time_i == time_j && status_i == status_j) {
        continue;
      }

      permissible += 1;

      if (time_i < time_j && sum_chf[i] > sum_chf[j]) {
        concordance += 1;
      } else if (time_j < time_i && sum_chf[j] > sum_chf[i]) {
        concordance += 1;
      } else if (sum_chf[i] == sum_chf[j]) {
        concordance += 0.5;
      }

    }
  }

  return (concordance / permissible);

}

std::string beautifyTime(uint seconds) {
  std::string result;

  // Add seconds, minutes, hours, days if larger than zero
  uint out_seconds = (uint) seconds % 60;
  result = std::to_string(out_seconds) + " seconds";
  uint out_minutes = (seconds / 60) % 60;
  if (seconds / 60 == 0) {
    return result;
  } else if (out_minutes == 1) {
    result = "1 minute, " + result;
  } else {
    result = std::to_string(out_minutes) + " minutes, " + result;
  }
  uint out_hours = (seconds / 3600) % 24;
  if (seconds / 3600 == 0) {
    return result;
  } else if (out_hours == 1) {
    result = "1 hour, " + result;
  } else {
    result = std::to_string(out_hours) + " hours, " + result;
  }
  uint out_days = (seconds / 86400);
  if (out_days == 0) {
    return result;
  } else if (out_days == 1) {
    result = "1 day, " + result;
  } else {
    result = std::to_string(out_days) + " days, " + result;
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

void shuffleAndSplit(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all, size_t n_first,
    std::mt19937_64 random_number_generator) {

  // Reserve space
  first_part.resize(n_all);

  // Fill with 0..n_all-1 and shuffle
  std::iota(first_part.begin(), first_part.end(), 0);
  std::shuffle(first_part.begin(), first_part.end(), random_number_generator);

  // Copy to second part
  second_part.resize(n_all - n_first);
  std::copy(first_part.begin() + n_first, first_part.end(), second_part.begin());

  // Resize first part
  first_part.resize(n_first);
}

