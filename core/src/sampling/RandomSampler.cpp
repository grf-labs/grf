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

#include "RandomSampler.h"

RandomSampler::RandomSampler(uint seed,
                                   SamplingOptions options):
    options(options) {
  random_number_generator.seed(seed);

}

void RandomSampler::sample(size_t num_samples,
                           double sample_fraction,
                           std::vector<size_t>& samples,
                           std::vector<size_t>& oob_samples) {
  bool sample_with_replacement = options.get_sample_with_replacement();
  if (options.get_case_weights().empty()) {
    if (sample_with_replacement) {
      bootstrap(num_samples, sample_fraction, samples, oob_samples);
    } else {
      bootstrap_without_replacement(num_samples, sample_fraction, samples, oob_samples);
    }
  } else {
    if (sample_with_replacement) {
      bootstrap_weighted(num_samples, sample_fraction, samples, oob_samples);
    } else {
      bootstrap_without_replacement_weighted(num_samples, sample_fraction, samples, oob_samples);
    }
  }
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

void RandomSampler::bootstrap(size_t num_samples,
                              double sample_fraction,
                              std::vector<size_t>& samples,
                              std::vector<size_t>& oob_sample) {
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;

  // Reserve space, reserve a little more to be safe
  samples.reserve(num_samples_inbag);
  oob_sample.reserve(num_samples * (exp(-sample_fraction) + 0.1));

  std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

  // Start with all samples OOB
  std::vector<size_t> inbag_counts;
  inbag_counts.resize(num_samples, 0);

  // Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = unif_dist(random_number_generator);
    samples.push_back(draw);
    ++inbag_counts[draw];
  }

// Save OOB samples
  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_sample.push_back(s);
    }
  }
}

void RandomSampler::bootstrap_weighted(size_t num_samples,
                                       double sample_fraction,
                                       std::vector<size_t>& samples,
                                       std::vector<size_t>& oob_samples) {
  const std::vector<double>& case_weights = options.get_case_weights();

  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;

// Reserve space, reserve a little more to be save)
  samples.reserve(num_samples_inbag);
  oob_samples.reserve(num_samples * (exp(-sample_fraction) + 0.1));

  std::discrete_distribution<> weighted_dist(case_weights.begin(), case_weights.end());

// Start with all samples OOB
  std::vector<size_t> inbag_counts;
  inbag_counts.resize(num_samples, 0);

// Draw num_samples samples with replacement (n out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = weighted_dist(random_number_generator);
    samples.push_back(draw);
    ++inbag_counts[draw];
  }

  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_samples.push_back(s);
    }
  }
}

void RandomSampler::bootstrap_without_replacement(size_t num_samples,
                                                  double sample_fraction,
                                                  std::vector<size_t>& samples,
                                                  std::vector<size_t>& oob_samples) {
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;
  shuffle_and_split(samples, oob_samples, num_samples, num_samples_inbag);
}

void RandomSampler::bootstrap_without_replacement_weighted(size_t num_samples,
                                                           double sample_fraction,
                                                           std::vector<size_t>& samples,
                                                           std::vector<size_t>& oob_samples) {
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;

  draw_without_replacement_weighted(samples,
                                    num_samples - 1,
                                    num_samples_inbag,
                                    options.get_case_weights());

  std::vector<size_t> inbag_counts;
  inbag_counts.resize(num_samples, 0);
  for (auto& sample : samples) {
    inbag_counts[sample] = 1;
  }

  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_samples.push_back(s);
    }
  }
}

void RandomSampler::shuffle_and_split(std::vector<size_t>& first_part,
                                      std::vector<size_t>& second_part,
                                      size_t n_all,
                                      size_t n_first) {
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
