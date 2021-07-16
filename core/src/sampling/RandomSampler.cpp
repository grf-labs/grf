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

#include <algorithm>
#include <random>

#include "RandomSampler.h"

namespace grf {

RandomSampler::RandomSampler(uint seed,
                             const SamplingOptions& options) :
    options(options) {
  random_number_generator.seed(seed);
}

void RandomSampler::sample_clusters(size_t num_rows,
                                    double sample_fraction,
                                    std::vector<size_t>& samples) {
  if (options.get_clusters().empty()) {
    sample(num_rows, sample_fraction, samples);
  } else {
    size_t num_samples = options.get_clusters().size();
    sample(num_samples, sample_fraction, samples);
  }
}

void RandomSampler::sample(size_t num_samples,
                           double sample_fraction,
                           std::vector<size_t>& samples) {
  size_t num_samples_inbag = static_cast<size_t>(num_samples * sample_fraction);
  shuffle_and_split(samples, num_samples, num_samples_inbag);
}

void RandomSampler::subsample(const std::vector<size_t>& samples,
                              double sample_fraction,
                              std::vector<size_t>& subsamples) {
  std::vector<size_t> shuffled_sample(samples);
  nonstd::shuffle(shuffled_sample.begin(), shuffled_sample.end(), random_number_generator);

  uint subsample_size = (uint) std::ceil(samples.size() * sample_fraction);
  subsamples.resize(subsample_size);
  std::copy(shuffled_sample.begin(),
            shuffled_sample.begin() + subsamples.size(),
            subsamples.begin());
}

void RandomSampler::subsample(const std::vector<size_t>& samples,
                              double sample_fraction,
                              std::vector<size_t>& subsamples,
                              std::vector<size_t>& oob_samples) {
  std::vector<size_t> shuffled_sample(samples);
  nonstd::shuffle(shuffled_sample.begin(), shuffled_sample.end(), random_number_generator);

  size_t subsample_size = (size_t) std::ceil(samples.size() * sample_fraction);
  subsamples.resize(subsample_size);
  oob_samples.resize(samples.size() - subsample_size);

  std::copy(shuffled_sample.begin(),
            shuffled_sample.begin() + subsamples.size(),
            subsamples.begin());
  std::copy(shuffled_sample.begin() + subsamples.size(),
            shuffled_sample.end(),
            oob_samples.begin());
}

void RandomSampler::subsample_with_size(const std::vector<size_t>& samples,
                                        size_t subsample_size,
                                        std::vector<size_t>& subsamples) {
  std::vector<size_t> shuffled_sample(samples);
  nonstd::shuffle(shuffled_sample.begin(), shuffled_sample.end(), random_number_generator);

  subsamples.resize(subsample_size);
  std::copy(shuffled_sample.begin(),
            shuffled_sample.begin() + subsamples.size(),
            subsamples.begin());
}

void RandomSampler::sample_from_clusters(const std::vector<size_t>& clusters,
                                         std::vector<size_t>& samples) {
  if (options.get_clusters().empty()) {
    samples = clusters;
  } else {
    const std::vector<std::vector<size_t>>& samples_by_cluster = options.get_clusters();
    for (size_t cluster : clusters) {
      const std::vector<size_t>& cluster_samples = samples_by_cluster[cluster];

      // Draw samples_per_cluster observations from each cluster. If the cluster is
      // smaller than the samples_per_cluster parameter, just use the whole cluster.
      if (cluster_samples.size() <= options.get_samples_per_cluster()) {
        samples.insert(samples.end(), cluster_samples.begin(), cluster_samples.end());
      } else {
        std::vector<size_t> subsamples;
        subsample_with_size(cluster_samples, options.get_samples_per_cluster(), subsamples);
        samples.insert(samples.end(), subsamples.begin(), subsamples.end());
      }
    }
  }
}

void RandomSampler::get_samples_in_clusters(const std::vector<size_t>& clusters,
                                            std::vector<size_t>& samples) {
  if (options.get_clusters().empty()) {
    samples = clusters;
  } else {
    for (size_t cluster : clusters) {
      const std::vector<size_t>& cluster_samples = options.get_clusters()[cluster];
      samples.insert(samples.end(), cluster_samples.begin(), cluster_samples.end());
    }
  }
}

void RandomSampler::shuffle_and_split(std::vector<size_t>& samples,
                                      size_t n_all,
                                      size_t size) {
  samples.resize(n_all);

  // Fill with 0..n_all-1 and shuffle
  std::iota(samples.begin(), samples.end(), 0);
  nonstd::shuffle(samples.begin(), samples.end(), random_number_generator);

  samples.resize(size);
}

void RandomSampler::draw(std::vector<size_t>& result,
                         size_t max,
                         const std::set<size_t>& skip,
                         size_t num_samples) {
  if (num_samples < max / 10) {
    draw_simple(result, max, skip, num_samples);
  } else {
    draw_fisher_yates(result, max, skip, num_samples);
  }
}

void RandomSampler::draw_simple(std::vector<size_t>& result,
                                size_t max,
                                const std::set<size_t>& skip,
                                size_t num_samples) {
  result.resize(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(max, false);

  nonstd::uniform_int_distribution<size_t> unif_dist(0, max - 1 - skip.size());
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
    result[i] = draw;
  }
}

void RandomSampler::draw_fisher_yates(std::vector<size_t>& result,
                                      size_t max,
                                      const std::set<size_t>& skip,
                                      size_t num_samples) {

  // Populate result vector with 0,...,max-1
  result.resize(max);
  std::iota(result.begin(), result.end(), 0);

  // Remove values that are to be skipped
  std::for_each(skip.rbegin(), skip.rend(),
                [&](size_t i) { result.erase(result.begin() + i); }
  );

  // Draw without replacement using Fisher Yates algorithm
  nonstd::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t i = 0; i < num_samples; ++i) {
    size_t j = static_cast<size_t>(i + distribution(random_number_generator) * (max - skip.size() - i));
    std::swap(result[i], result[j]);
  }

  result.resize(num_samples);
}

size_t RandomSampler::sample_poisson(size_t mean) {
  nonstd::poisson_distribution<size_t> distribution(static_cast<double>(mean));
  return distribution(random_number_generator);
}

} // namespace grf
