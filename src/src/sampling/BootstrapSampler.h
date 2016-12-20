
#ifndef RANGER_BOOTSTRAPSAMPLER_H
#define RANGER_BOOTSTRAPSAMPLER_H

#include <cstddef>
#include <vector>
#include <random>
#include "globals.h"

class BootstrapSampler {
public:
  BootstrapSampler(size_t num_samples,
                   uint seed,
                   bool sample_with_replacement,
                   double sample_fraction,
                   bool keep_inbag,
                   std::vector<double> *case_weights);

  void sample(std::vector<std::vector<size_t>>& sampleIDs);

  const std::vector<size_t>& getOobSampleIDs() const {
    return oob_sampleIDs;
  }

  const size_t getNumSamplesOob() const {
    return num_samples_oob;
  };

  const std::vector<size_t>& getInbagCounts() const {
    return inbag_counts;
  }

  /**
   * Draw random numbers in a range without replacement and skip values.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void drawWithoutReplacementSkip(std::vector<size_t> &result,
                                  size_t range_length,
                                  std::vector<size_t> &skip,
                                  size_t num_samples);

  /**
   * Draw random numers without replacement and with weighted probabilites from vector of indices.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param indices Vector with numbers to draw
   * @param num_samples Number of samples to draw
   * @param weights A weight for each element of indices
   */
  void drawWithoutReplacementWeighted(std::vector<size_t> &result,
                                      std::vector<size_t> &indices,
                                      size_t num_samples,
                                      std::vector<double> &weights);

  /**
   * Draw random numers without replacement and with weighted probabilites from vector of indices.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param random_number_generator Random number generator
   * @param indices Vector with numbers to draw
   * @param num_samples Number of samples to draw
   * @param weights A weight for each element of indices
   */
  void drawWithoutReplacementWeighted(std::vector<size_t> &result,
                                      size_t max_index,
                                      size_t num_samples,
                                      std::vector<double> &weights);

  /**
   * Create numbers from 0 to n_all-1, shuffle and split in two parts.
   * @param first_part First part
   * @param second_part Second part
   * @param n_all Number elements
   * @param n_first Number of elements of first part
   */
  void shuffleAndSplit(std::vector<size_t> &first_part,
                       std::vector<size_t> &second_part,
                       size_t n_all,
                       size_t n_first);

private:
  void bootstrap(std::vector<std::vector<size_t>>& sampleIDs);
  void bootstrapWithoutReplacement(std::vector<std::vector<size_t>>& sampleIDs);
  void bootstrapWeighted(std::vector<std::vector<size_t>>& sampleIDs);
  void bootstrapWithoutReplacementWeighted(std::vector<std::vector<size_t>>& sampleIDs);

  /**
   * Simple algorithm for sampling without replacement, faster for smaller num_samples
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void drawWithoutReplacementSimple(std::vector<size_t> &result,
                                    size_t max,
                                    std::vector<size_t> &skip,
                                    size_t num_samples);

  /**
   * Knuth's algorithm for sampling without replacement, faster for larger num_samples
   * Idea from Knuth 1985, The Art of Computer Programming, Vol. 2, Sec. 3.4.2 Algorithm S
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void drawWithoutReplacementKnuth(std::vector<size_t> &result,
                                   size_t max,
                                   std::vector<size_t> &skip,
                                   size_t num_samples);

  bool sample_with_replacement;
  double sample_fraction;
  bool keep_inbag;
  std::vector<size_t> inbag_counts;
  std::vector<double>* case_weights;
  std::vector<size_t> oob_sampleIDs;
  size_t num_samples_oob;
  std::mt19937_64 random_number_generator;

  size_t num_samples;
};


#endif //RANGER_BOOTSTRAPSAMPLER_H
