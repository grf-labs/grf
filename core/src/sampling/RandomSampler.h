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
                const SamplingOptions& options);

  /**
   * Samples some number of clusters, given the configuration in {@link SampleOptions}.
   *
   * If sample clustering is enabled, this method will return cluster IDs. Otherwise,
   * the returned IDs represent individual sample IDs for efficiency. They should still
   * be thought of as degenerate 'cluster IDs', and even if clustering is not enabled,
   * these IDs can be passed to the other cluster methods below.
   *
   * @param num_rows The total number of rows in the input data.
   * @param sample_fraction The fraction of clusters that should be in the sample.
   * @param samples An empty vector, which this method will populate with the sampled cluster IDs.
   */
  void sample_clusters(size_t num_rows,
                       double sample_fraction,
                       std::vector<size_t>& samples);

  /**
   * If clustering is enabled, draws the appropriate number of samples from the provided
   * cluster IDs. Otherwise, it is a no-op: we simply return the passed in 'cluster IDs',
   * as they already represent individual sample IDs.
   *
   * Note that the number of samples drawn is configured through
   * {@link SamplingOptions#get_samples_per_cluster}.
   */
  void sample_from_clusters(const std::vector<size_t>& clusters,
                            std::vector<size_t>& samples);

  /**
   * If clustering is enabled, returns all samples in the givenclusters. Otherwise,
   * Returns IDs from the provided cluster IDs. Otherwise, we return the passed-in
   * 'cluster IDs', as they already represent individual sample IDs.
   */
  void get_samples_in_clusters(const std::vector<size_t>& clusters,
                               std::vector<size_t>& samples);

  void sample(size_t num_samples,
              double sample_fraction,
              std::vector<size_t>& samples);

  void subsample(const std::vector<size_t>& samples,
                 double sample_fraction,
                 std::vector<size_t>& subsamples);

  void subsample(const std::vector<size_t>& samples,
                 double sample_fraction,
                 std::vector<size_t>& subsamples,
                 std::vector<size_t>& oob_samples);

  void subsample_with_size(const std::vector<size_t>& samples,
                           size_t subsample_size,
                           std::vector<size_t>& subsamples);

  /**
   * Draw random numbers in a range without replacement and skip values.
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param max Specifies the interval to draw from:  0 ... (max-1).
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void draw(std::vector<size_t>& result,
            size_t max,
            const std::set<size_t>& skip,
            size_t num_samples);

  size_t sample_poisson(size_t mean);

private:
 /**
  * Create numbers from 0 to n_all-1, then shuffle and select the first 'size' elements.
  *
  * @param samples A list of first 'size'n_first shuffled numbers
  * @param n_all Number elements
  * @param size Number of elements of to select
  */
  void shuffle_and_split(std::vector<size_t>& samples,
                         size_t n_all,
                         size_t size);

  /**
   * Simple algorithm for sampling without replacement, faster for smaller num_samples
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param range_length Length of range. Interval to draw from: 0..max-1
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void draw_simple(std::vector<size_t>& result,
                   size_t max,
                   const std::set<size_t>& skip,
                   size_t num_samples);

  /**
  * Draw random numbers without replacement and with weighted probabilites from vector of indices.
  * @param result Vector to add results to. Will not be cleaned before filling.
  * @param max Specifies the interval to draw from:  0 ... (max-1).
  * @param num_samples Number of samples to draw
  * @param weights A weight for each element of indices
  */
  void draw_weighted(std::vector<size_t>& result,
                     size_t max,
                     size_t num_samples,
                     const std::vector<double>& weights);

  /**
   * Knuth's algorithm for sampling without replacement, faster for larger num_samples
   * Idea from Knuth 1985, The Art of Computer Programming, Vol. 2, Sec. 3.4.2 Algorithm S
   * @param result Vector to add results to. Will not be cleaned before filling.
   * @param max Specifies the interval to draw from:  0 ... (max-1).
   * @param skip Values to skip
   * @param num_samples Number of samples to draw
   */
  void draw_knuth(std::vector<size_t>& result,
                  size_t max,
                  const std::set<size_t>& skip,
                  size_t num_samples);

  SamplingOptions options;
  std::mt19937_64 random_number_generator;
};


#endif //GRF_BOOTSTRAPSAMPLER_H
