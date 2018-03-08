/*-------------------------------------------------------------------------------
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

#ifndef GRF_BOOTSTRAPSAMPLER_H
#define GRF_BOOTSTRAPSAMPLER_H


#include "commons/globals.h"
#include "commons/utility.h"
#include "SamplingOptions.h"

#include <cstddef>
#include <random>
#include <set>
#include <vector>
#include <unordered_map>

class RandomSampler {
public:
  RandomSampler(uint seed,
                SamplingOptions options);

  void sample(size_t num_samples,
              double sample_fraction,
              std::vector<size_t>& samples,
              std::vector<size_t>& oob_samples);

  void sample_for_ci(Data* data,
                     double sample_fraction,
                     std::vector<size_t>& samples,
                     std::vector<size_t>& oob_samples);

  void subsample(const std::vector<size_t>& samples,
                 double sample_fraction,
                 std::vector<size_t>& subsamples,
                 std::vector<size_t>& oob_samples);

  /**
   * Split sample into subsample and oob_sample sets. The split is performed by cluster so that all observations
   * corresponding to an individual cluster are entirely contained in one of the two subsets.
   *
   * @param samples Observations to sample
   * @param sample_fraction Fraction of clusters that go into subsamples
   * @param subsamples Empty vector to be populated by subsample
   * @param oob_samples Empty vector to be populated by oob sample
   * @param clusters Cluster information for each entry in samples
   */
  void subsample_with_clusters(const std::vector<size_t>& samples,
                               double sample_fraction,
                               std::vector<size_t>& subsamples,
                               std::vector<size_t>& oob_samples,
                               std::vector<uint>& clusters);

  /**
   * Draw random numbers in a range without replacement and skip values.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void draw_without_replacement_skip(std::vector<size_t>& result,
                                     size_t range_length,
                                     const std::set<size_t>& skip,
                                     size_t num_samples);

  /**
   * Draw random numebrs without replacement and with weighted probabilites from vector of indices.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param indices Vector with numbers to draw
   * @param num_samples Number of samples to draw
   * @param weights A weight for each element of indices
   */
  void draw_without_replacement_weighted(std::vector<size_t>& result,
                                         const std::vector<size_t>& indices,
                                         size_t num_samples,
                                         const std::vector<double>& weights);

  /**
   * Draw random numbers without replacement and with weighted probabilites from vector of indices.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param random_number_generator Random number generator
   * @param indices Vector with numbers to draw
   * @param num_samples Number of samples to draw
   * @param weights A weight for each element of indices
   */
  void draw_without_replacement_weighted(std::vector<size_t>& result,
                                         size_t max_index,
                                         size_t num_samples,
                                         const std::vector<double>& weights);

  /**
   * Create numbers from 0 to n_all-1, shuffle and split in two parts.
   * @param first_part First part
   * @param second_part Second part
   * @param n_all Number elements
   * @param n_first Number of elements of first part
   */
  void shuffle_and_split(std::vector<size_t>& first_part,
                         std::vector<size_t>& second_part,
                         size_t n_all,
                         size_t n_first);

  size_t sample_poisson(size_t mean);

  /**
   * Performs bootstrapping at the cluster level rather than at the individual level.
   * @param num_samples Number of rows to be sampled from
   * @param sample_fraction Fraction to include inbag
   * @param samples Samples to be populated
   * @param oob_samples OOB samples to be populated
   * @param clusters Cluster information for each observation
   */
  void bootstrap_with_clusters(size_t num_samples,
                               double sample_fraction,
                               std::vector<size_t>& samples,
                               std::vector<size_t>& oob_samples,
                               std::unordered_map<uint, std::vector<size_t>>& clusters);

private:
  void bootstrap(size_t num_samples,
                 double sample_fraction,
                 std::vector<size_t>& samples,
                 std::vector<size_t>& oob_sample);

  void bootstrap_without_replacement(size_t num_samples,
                                     double sample_fraction,
                                     std::vector<size_t>& samples,
                                     std::vector<size_t>& oob_samples);

  void bootstrap_weighted(size_t num_samples,
                          double sample_fraction,
                          std::vector<size_t>& samples,
                          std::vector<size_t>& oob_samples);

  void bootstrap_without_replacement_weighted(size_t num_samples,
                                              double sample_fraction,
                                              std::vector<size_t>& samples,
                                              std::vector<size_t>& oob_samples);

  void bootstrap_without_oob(size_t num_samples,
                             double sample_fraction,
                             std::vector<size_t>& samples);

  /**
   * Simple algorithm for sampling without replacement, faster for smaller num_samples
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void draw_without_replacement(std::vector<size_t>& result,
                                size_t max,
                                const std::set<size_t>& skip,
                                size_t num_samples);

  /**
   * Knuth's algorithm for sampling without replacement, faster for larger num_samples
   * Idea from Knuth 1985, The Art of Computer Programming, Vol. 2, Sec. 3.4.2 Algorithm S
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void draw_without_replacement_knuth(std::vector<size_t>& result,
                                      size_t max,
                                      const std::set<size_t>& skip,
                                      size_t num_samples);

  SamplingOptions options;
  std::mt19937_64 random_number_generator;
};


#endif //GRF_BOOTSTRAPSAMPLER_H
