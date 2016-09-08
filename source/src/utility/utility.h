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

#ifdef R_BUILD
#include <Rinternals.h>
#endif

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
 * Draw random numbers in a range without replacement and skip values.
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param skip Values to skip
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementSkip(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t range_length, std::vector<size_t>& skip, size_t num_samples);

/**
 * Simple algorithm for sampling without replacement, faster for smaller num_samples
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param skip Values to skip
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    std::vector<size_t>& skip, size_t num_samples);

/**
 * Knuth's algorithm for sampling without replacement, faster for larger num_samples
 * Idea from Knuth 1985, The Art of Computer Programming, Vol. 2, Sec. 3.4.2 Algorithm S
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param skip Values to skip
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementKnuth(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    std::vector<size_t>& skip, size_t num_samples);

/**
 * Draw random numers without replacement and with weighted probabilites from vector of indices.
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param indices Vector with numbers to draw
 * @param num_samples Number of samples to draw
 * @param weights A weight for each element of indices
 */
void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    std::vector<size_t>& indices, size_t num_samples, std::vector<double>& weights);

/**
 * Draw random numers without replacement and with weighted probabilites from 0..n-1.
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param max_index Maximum index to draw
 * @param num_samples Number of samples to draw
 * @param weights A weight for each element of indices
 */
void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max_index, size_t num_samples, std::vector<double>& weights);

/**
 * Returns the most frequent class index of a vector with counts for the classes. Returns a random class if counts are equal.
 * @param class_count Vector with class counts
 * @param random_number_generator Random number generator
 * @return Most frequent class index. Out of range index if all 0.
 */
template<typename T>
size_t mostFrequentClass(std::vector<T>& class_count, std::mt19937_64 random_number_generator) {
  std::vector<size_t> major_classes;

  // Find maximum count
  T max_count = 0;
  for (size_t i = 0; i < class_count.size(); ++i) {
    T count = class_count[i];
    if (count > max_count) {
      max_count = count;
      major_classes.clear();
      major_classes.push_back(i);
    } else if (count == max_count) {
      major_classes.push_back(i);
    }
  }

  if (max_count == 0) {
    return class_count.size();
  } else if (major_classes.size() == 1) {
    return major_classes[0];
  } else {
    // Choose randomly
    std::uniform_int_distribution<size_t> unif_dist(0, major_classes.size() - 1);
    return major_classes[unif_dist(random_number_generator)];
  }
}

/**
 * Returns the most frequent value of a map with counts for the values. Returns a random class if counts are equal.
 * @param class_count Map with classes and counts
 * @param random_number_generator Random number generator
 * @return Most frequent value
 */
double mostFrequentValue(std::unordered_map<double, size_t>& class_count, std::mt19937_64 random_number_generator);

/**
 * Compute concordance index for given data and summed cumulative hazard function/estimate
 * @param data Pointer to Data object
 * @param sum_chf Summed chf over timepoints for each sample
 * @param dependent_varID ID of dependent variable
 * @param status_varID ID of status variable
 * @param sample_IDs IDs of samples, for example OOB samples
 * @return concordance index
 */
double computeConcordanceIndex(Data* data, std::vector<double>& sum_chf, size_t dependent_varID, size_t status_varID,
    std::vector<size_t>& sample_IDs);

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
 * Create numbers from 0 to n_all-1, shuffle and split in two parts.
 * @param first_part First part
 * @param second_part Second part
 * @param n_all Number elements
 * @param n_first Number of elements of first part
 * @param random_number_generator Random number generator
 */
void shuffleAndSplit(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all, size_t n_first,
    std::mt19937_64 random_number_generator);

/**
 * Check if not too many factor levels and all values in unordered categorical variables are positive integers.
 * @param data Pointer to data object
 * @param unordered_variable_names Names of unordered variables
 * @return Error message, empty if no problem occured
 */
std::string checkUnorderedVariables(Data* data, std::vector<std::string> unordered_variable_names);

/**
 * Check if all values in double vector are positive integers.
 * @param all_values Double vector to check
 * @return True if all values are positive integers
 */
bool checkPositiveIntegers(std::vector<double>& all_values);

/**
 * Compute p-value for maximally selected rank statistics using Lau92 approximation
 * See Lausen, B. & Schumacher, M. (1992). Biometrics 48, 73-85.
 * @param b Quantile
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 * @return p-value for quantile b
 */
double maxstatPValueLau92(double b, double minprop, double maxprop);

/**
 * Compute p-value for maximally selected rank statistics using Lau94 approximation
 * See Lausen, B., Sauerbrei, W. & Schumacher, M. (1994). Computational Statistics. 483-496.
 * @param b Quantile
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 * @param N Number of observations
 * @param m Vector with number of observations smaller or equal than cutpoint, sorted, only for unique cutpoints
 * @return p-value for quantile b
 */
double maxstatPValueLau94(double b, double minprop, double maxprop, size_t N, std::vector<size_t>& m);

/**
 * Standard normal density
 * @param x Quantile
 * @return Standard normal density at quantile x
 */
double dstdnorm(double x);

/**
 * Standard normal distribution
 * @param x Quantile
 * @return Standard normal distribution at quantile x
 */
double pstdnorm(double x);

/**
 * Adjust p-values with Benjamini/Hochberg
 * @param unadjusted_pvalues Unadjusted p-values (input)
 * @param adjusted_pvalues Adjusted p-values (result)
 */
std::vector<double> adjustPvalues(std::vector<double>& unadjusted_pvalues);

/**
 * Get indices of sorted values
 * @param values Values to sort
 * @param decreasing Order decreasing
 * @return Indices of sorted values
 */
template<typename T>
std::vector<size_t> order(std::vector<T>& values, bool decreasing) {
  // Create index vector
  std::vector<size_t> indices(values.size());
  std::iota(indices.begin(), indices.end(), 0);

  // Sort index vector based on value vector
  if (decreasing) {
    std::sort(std::begin(indices), std::end(indices), [&](size_t i1, size_t i2) {return values[i1] > values[i2];});
  } else {
    std::sort(std::begin(indices), std::end(indices), [&](size_t i1, size_t i2) {return values[i1] < values[i2];});
  }
  return indices;
}

/**
 * Get indices of sorted values, compute in-place of Data.
 * @param data Data object
 * @param sampleIDs IDs of samples to sort
 * @param varID ID of variable to sort by
 * @param decreasing Order decreasing
 * @return Indices of sorted values
 */
std::vector<size_t> orderInData(Data* data, std::vector<size_t>& sampleIDs, size_t varID, bool decreasing);

/**
 * Sample ranks starting from 1. Ties are given the average rank.
 * @param values Values to rank
 * @return Ranks of input values
 */
template<typename T>
std::vector<double> rank(std::vector<T>& values) {
  size_t num_values = values.size();

  // Order
  std::vector<size_t> indices = order(values, false);

  // Compute ranks, start at 1
  std::vector<double> ranks(num_values);
  size_t reps = 1;
  for (size_t i = 0; i < num_values; i += reps) {

    // Find number of replications
    reps = 1;
    while (i + reps < num_values && values[indices[i]] == values[indices[i + reps]]) {
      ++reps;
    }

    // Assign rank to all replications
    for (size_t j = 0; j < reps; ++j)
      ranks[indices[i + j]] = (2 * (double) i + (double) reps - 1) / 2 + 1;
  }

  return ranks;
}

/**
 * Compute Logrank scores for survival times
 * @param time Survival time
 * @param status Censoring indicator
 * @return Logrank scores
 */
std::vector<double> logrankScores(std::vector<double>& time, std::vector<double>& status);

/**
 * Compute Logrank scores for survival times directly from Data object
 * @param data Pointer to Data object
 * @param time_varID variable ID for time column
 * @param status_varID Variable ID for status column
 * @param sampleIDs IDs of samples to include
 * @return Logrank scores
 */
std::vector<double> logrankScoresData(Data* data, size_t time_varID, size_t status_varID,
    std::vector<size_t> sampleIDs);

/**
 * Compute maximally selected rank statistics
 * @param scores Scores for dependent variable (y)
 * @param x Independent variable
 * @param indices Ordering of x values
 * @param best_maxstat Maximally selected statistic (output)
 * @param best_split_value Split value for maximally selected statistic (output)
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 */
void maxstat(std::vector<double>& scores, std::vector<double>& x, std::vector<size_t>& indices, double& best_maxstat,
    double& best_split_value, double minprop, double maxprop);

/**
 * Compute maximally selected rank statistics in Data object
 * @param scores Scores for dependent variable (y)
 * @param data Data object
 * @param sampleIDs IDs of samples to consider
 * @param varID ID of variable to consider
 * @param indices Ordering of x values
 * @param best_maxstat Maximally selected statistic (output)
 * @param best_split_value Split value for maximally selected statistic (output)
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 */
void maxstatInData(std::vector<double>& scores, Data* data, std::vector<size_t>& sampleIDs, size_t varID,
    std::vector<size_t>& indices, double& best_maxstat, double& best_split_value, double minprop, double maxprop);

/**
 * Compute number of samples smaller or equal than each unique value in x
 * @param x Value vector
 * @param indices Ordering of x
 * @return Vector of number of samples smaller or equal than each unique value in x
 */
std::vector<size_t> numSamplesLeftOfCutpoint(std::vector<double>& x, std::vector<size_t>& indices);

/**
 * Compute number of samples smaller or equal than each unique value in data
 * @param data Data object
 * @param sampleIDs IDs of samples to consider
 * @param varID ID of variable to consider
 * @param indices Ordering of x
 * @return Vector of number of samples smaller or equal than each unique value in x
 */
std::vector<size_t> numSamplesLeftOfCutpointInData(Data* data, std::vector<size_t>& sampleIDs, size_t varID,
    std::vector<size_t>& indices);

// User interrupt from R
#ifdef R_BUILD
static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

inline bool checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}
#endif

#endif /* UTILITY_H_ */
