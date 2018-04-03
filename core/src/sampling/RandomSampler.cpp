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
#include <unordered_map>

#include "RandomSampler.h"

RandomSampler::RandomSampler(uint seed,
                             SamplingOptions options):
    options(options) {
  random_number_generator.seed(seed);

}

void RandomSampler::sample(size_t num_samples,
                           double sample_fraction,
                           std::vector<size_t>& samples) {
  bool sample_with_replacement = options.get_sample_with_replacement();
  if (options.get_sample_weights().empty()) {
    if (sample_with_replacement) {
      bootstrap(num_samples, sample_fraction, samples);
    } else {
      bootstrap_without_replacement(num_samples, sample_fraction, samples);
    }
  } else {
    if (sample_with_replacement) {
      bootstrap_weighted(num_samples, sample_fraction, samples);
    } else {
      bootstrap_without_replacement_weighted(num_samples, sample_fraction, samples);
    }
  }
}

void RandomSampler::sample_for_ci(Data* data,
                                  double sample_fraction,
                                  std::vector<size_t>& samples,
                                  std::vector<size_t>& oob_samples) {
  size_t num_samples = data->get_num_rows();
  uint samples_per_cluster = options.get_samples_per_cluster();
  if (samples_per_cluster > 1) {
    bootstrap_with_clusters(num_samples, sample_fraction, samples, oob_samples, data->get_cluster_map());
  } else {
    sample(num_samples, sample_fraction, samples);
  }
}

void RandomSampler::subsample_for_ci(const std::vector<size_t>& samples,
                                     double sample_fraction,
                                     std::vector<size_t>& subsamples,
                                     std::vector<size_t>& oob_samples,
                                     Data* data) {
  if (options.get_samples_per_cluster() > 1) {
    auto clusters = data->get_clusters();
    subsample_with_clusters(samples, sample_fraction, subsamples, oob_samples, clusters);
  } else {
    subsample(samples, sample_fraction, subsamples, oob_samples);
  }
}

void RandomSampler::subsample(const std::vector<size_t>& samples,
                              double sample_fraction,
                              std::vector<size_t>& subsamples) {
  std::vector<size_t> shuffled_sample(samples);
  std::shuffle(shuffled_sample.begin(), shuffled_sample.end(), random_number_generator);

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
  std::shuffle(shuffled_sample.begin(), shuffled_sample.end(), random_number_generator);

  uint subsample_size = (uint) std::ceil(samples.size() * sample_fraction);
  subsamples.resize(subsample_size);
  oob_samples.resize(samples.size() - subsample_size);

  std::copy(shuffled_sample.begin(),
            shuffled_sample.begin() + subsamples.size(),
            subsamples.begin());
  std::copy(shuffled_sample.begin() + subsamples.size(),
            shuffled_sample.end(),
            oob_samples.begin());
}

void RandomSampler::subsample_with_clusters(const std::vector<size_t>& samples,
                                            double sample_fraction,
                                            std::vector<size_t>& subsamples,
                                            std::vector<size_t>& oob_samples,
                                            std::vector<uint>& clusters) {
  // get map of cluster to observations
  std::unordered_map<uint, std::vector<size_t>> cluster_map;
  for (auto const& obs : samples) {
    auto cluster = clusters[obs];
    cluster_map[cluster].push_back(obs);
  }

  // Now shuffle the clusters
  std::vector<uint> shuffled_clusters;
  shuffled_clusters.reserve(cluster_map.size());
  for(auto const& imap: cluster_map) {
    shuffled_clusters.push_back(imap.first);
  }

  std::shuffle(shuffled_clusters.begin(), shuffled_clusters.end(), random_number_generator);

  // finally populate vectors with proper observations
  uint subsample_size = (uint) std::ceil(shuffled_clusters.size() * sample_fraction);

  for (size_t s = 0; s < shuffled_clusters.size(); ++s) {
    auto cluster = shuffled_clusters[s];
    if (s < subsample_size) {
      subsamples.insert(subsamples.end(), cluster_map[cluster].begin(), cluster_map[cluster].end());
    } else {
      oob_samples.insert(oob_samples.end(), cluster_map[cluster].begin(), cluster_map[cluster].end());
    }
  }
}

void RandomSampler::bootstrap(size_t num_samples,
                              double sample_fraction,
                              std::vector<size_t>& samples) {

  // Reserve space, reserve a little more to be safe
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;
  samples.reserve(num_samples_inbag);

  std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

  // Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = unif_dist(random_number_generator);
    samples.push_back(draw);
  }
}

void RandomSampler::bootstrap_weighted(size_t num_samples,
                                       double sample_fraction,
                                       std::vector<size_t>& samples) {
  const std::vector<double>& sample_weights = options.get_sample_weights();

  // Reserve space, reserve a little more to be save)
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;
  samples.reserve(num_samples_inbag);

  std::discrete_distribution<> weighted_dist(sample_weights.begin(), sample_weights.end());

  // Draw num_samples samples with replacement (n out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = weighted_dist(random_number_generator);
    samples.push_back(draw);
  }
}

void RandomSampler::bootstrap_without_replacement(size_t num_samples,
                                                  double sample_fraction,
                                                  std::vector<size_t>& samples) {
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;
  shuffle_and_split(samples, num_samples, num_samples_inbag);
}

void RandomSampler::bootstrap_without_replacement_weighted(size_t num_samples,
                                                           double sample_fraction,
                                                           std::vector<size_t>& samples) {
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;
  draw_without_replacement_weighted(samples,
                                    num_samples - 1,
                                    num_samples_inbag,
                                    options.get_sample_weights());
}

void RandomSampler::bootstrap_without_oob(size_t num_samples,
                                          double sample_fraction,
                                          std::vector<size_t>& samples) {
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;

  // Reserve space, reserve a little more to be safe
  samples.reserve(num_samples_inbag);

  std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

  // Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = unif_dist(random_number_generator);
    samples.push_back(draw);
  }
}

void RandomSampler::bootstrap_with_clusters(size_t num_samples,
                                            double sample_fraction,
                                            std::vector<size_t>& samples,
                                            std::vector<size_t>& oob_sample,
                                            std::unordered_map<uint, std::vector<size_t>>& clusters) {

  // Reserve space
  size_t num_clusters = clusters.size();
  size_t num_samples_inbag = (size_t) num_clusters * sample_fraction * options.get_samples_per_cluster();
  samples.reserve(num_samples_inbag);
  oob_sample.reserve(num_samples * (std::exp(-sample_fraction) + 0.1));

  // First sample at the cluster level
  std::vector<size_t> cluster_samples;
  RandomSampler::bootstrap_without_oob(num_clusters, sample_fraction, cluster_samples);

  // Count number of times each cluster is sampled
  std::unordered_map<size_t, size_t> inbag_cluster_counts;
  for (size_t sampled_cluster : cluster_samples) {
    inbag_cluster_counts[sampled_cluster]++;
  }

  // Add all non-sampled clusters to oob
  for (size_t i = 0; i < num_clusters; ++i) {
    if (inbag_cluster_counts.find(i) == inbag_cluster_counts.end()) {
      oob_sample.insert(oob_sample.end(), clusters[i].begin(), clusters[i].end());
    }
  }

  // Sample from each selected cluster
  for (auto const& entry : inbag_cluster_counts) {
    std::vector<size_t> cluster_obs = clusters[entry.first];
    size_t cluster_size = cluster_obs.size();

    // Case where same cluster sampled multiple times
    for (size_t i = 0; i < entry.second; ++i) {
      std::vector<size_t> sample_indices;
      double within_cluster_sample_fraction = (double) options.get_samples_per_cluster() / cluster_size;
      RandomSampler::bootstrap_without_oob(cluster_size,
                                           within_cluster_sample_fraction,
                                           sample_indices);

      for (size_t k : sample_indices) {
        samples.push_back(cluster_obs[k]);
      }
    }
  }
}

void RandomSampler::shuffle_and_split(std::vector<size_t> &samples,
                                      size_t n_all,
                                      size_t size) {
  samples.resize(n_all);

  // Fill with 0..n_all-1 and shuffle
  std::iota(samples.begin(), samples.end(), 0);
  std::shuffle(samples.begin(), samples.end(), random_number_generator);

  samples.resize(size);
}

void RandomSampler::draw_without_replacement_skip(std::vector<size_t>& result,
                                                  size_t max,
                                                  const std::set<size_t>& skip,
                                                  size_t num_samples) {
  if (num_samples < max / 2) {
    draw_without_replacement(result, max, skip, num_samples);
  } else {
    draw_without_replacement_knuth(result, max, skip, num_samples);
  }
}

void RandomSampler::draw_without_replacement(std::vector<size_t>& result,
                                             size_t max,
                                             const std::set<size_t>& skip,
                                             size_t num_samples) {
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

void RandomSampler::draw_without_replacement_knuth(std::vector<size_t>& result,
                                                   size_t max,
                                                   const std::set<size_t>& skip,
                                                   size_t num_samples) {
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

void RandomSampler::draw_without_replacement_weighted(std::vector<size_t>& result,
                                                      const std::vector<size_t>& indices,
                                                      size_t num_samples,
                                                      const std::vector<double>& weights) {
  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(indices.size(), false);

  std::discrete_distribution<> weighted_dist(weights.begin(), weights.end());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = weighted_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(indices[draw]);
  }
}

void RandomSampler::draw_without_replacement_weighted(std::vector<size_t>& result,
                                                      size_t max_index,
                                                      size_t num_samples,
                                                      const std::vector<double>& weights) {
  result.reserve(num_samples);

  // Set all to not selected
  std::vector<bool> temp;
  temp.resize(max_index + 1, false);

  std::discrete_distribution<> weighted_dist(weights.begin(), weights.end());
  for (size_t i = 0; i < num_samples; ++i) {
    size_t draw;
    do {
      draw = weighted_dist(random_number_generator);
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

size_t RandomSampler::sample_poisson(size_t mean) {
  std::poisson_distribution<size_t> distribution(mean);
  return distribution(random_number_generator);
}
