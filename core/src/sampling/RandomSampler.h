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
              std::vector<size_t>& samples);

  void subsample(const std::vector<size_t>& samples,
                 double sample_fraction,
                 std::vector<size_t>& subsamples);

  /**
   * Sampling to be used only in the case of calculating confidence intervals. If clustered standard errors
   * are enabled, sampling will return the ids corresponding to sample clusters. Otherwise, ids will correspond
   * to individual observations.
   *
   * @param data
   * @param sample_fraction
   * @param samples
   */
  void sample_for_ci(Data* data,
                     double sample_fraction,
                     std::vector<size_t>& samples);

  void subsample(const std::vector<size_t>& samples,
                 double sample_fraction,
                 std::vector<size_t>& subsamples,
                 std::vector<size_t>& oob_samples);

  void sample_from_clusters(const std::vector<size_t>& cluster_samples,
                            std::vector<size_t>& samples);

  void get_oob_from_clusters(const std::vector<size_t>& cluster_oob_sample,
                             std::vector<size_t>& oob_samples);

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

  bool clustering_enabled() const;

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
