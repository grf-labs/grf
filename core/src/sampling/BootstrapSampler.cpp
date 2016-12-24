
#include "BootstrapSampler.h"

BootstrapSampler::BootstrapSampler(size_t num_samples,
                                   uint seed,
                                   SamplingOptions* options):
    num_samples(num_samples),
    options(options) {
  random_number_generator.seed(seed);
}

void BootstrapSampler::sample(std::vector<std::vector<size_t>>& sampleIDs) {
  bool sample_with_replacement = options->get_sample_with_replacement();
  if (options->get_case_weights()->empty()) {
    if (sample_with_replacement) {
      bootstrap(sampleIDs);
    } else {
      bootstrapWithoutReplacement(sampleIDs);
    }
  } else {
    if (sample_with_replacement) {
      bootstrapWeighted(sampleIDs);
    } else {
      bootstrapWithoutReplacementWeighted(sampleIDs);
    }
  }
}

void BootstrapSampler::bootstrap(std::vector<std::vector<size_t>>& sampleIDs) {
  double sample_fraction = options->get_sample_fraction();
  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;

// Reserve space, reserve a little more to be save)
  sampleIDs[0].reserve(num_samples_inbag);
  oob_sampleIDs.reserve(num_samples * (exp(-sample_fraction) + 0.1));

  std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

// Start with all samples OOB
  inbag_counts.resize(num_samples, 0);

// Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = unif_dist(random_number_generator);
    sampleIDs[0].push_back(draw);
    ++inbag_counts[draw];
  }

// Save OOB samples
  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_sampleIDs.push_back(s);
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  if (!options->get_keep_in_bag()) {
    inbag_counts.clear();
  }
}

void BootstrapSampler::bootstrapWeighted(std::vector<std::vector<size_t>>& sampleIDs) {
  double sample_fraction = options->get_sample_fraction();
  std::vector<double>* case_weights = options->get_case_weights();

  size_t num_samples_inbag = (size_t) num_samples * sample_fraction;

// Reserve space, reserve a little more to be save)
  sampleIDs[0].reserve(num_samples_inbag);
  oob_sampleIDs.reserve(num_samples * (exp(-sample_fraction) + 0.1));

  std::discrete_distribution<> weighted_dist(case_weights->begin(), case_weights->end());

// Start with all samples OOB
  inbag_counts.resize(num_samples, 0);

// Draw num_samples samples with replacement (n out of n) as inbag and mark as not OOB
  for (size_t s = 0; s < num_samples_inbag; ++s) {
    size_t draw = weighted_dist(random_number_generator);
    sampleIDs[0].push_back(draw);
    ++inbag_counts[draw];
  }

  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_sampleIDs.push_back(s);
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  if (!options->get_keep_in_bag()) {
    inbag_counts.clear();
  }
}

void BootstrapSampler::bootstrapWithoutReplacement(std::vector<std::vector<size_t>> &sampleIDs) {
  // Use fraction (default 63.21%) of the samples
  size_t num_samples_inbag = (size_t) num_samples * options->get_sample_fraction();
  shuffleAndSplit(sampleIDs[0], oob_sampleIDs, num_samples, num_samples_inbag);
  num_samples_oob = oob_sampleIDs.size();

  if (options->get_keep_in_bag()) {
    // All observation are 0 or 1 times inbag
    inbag_counts.resize(num_samples, 1);
    for (size_t i = 0; i < oob_sampleIDs.size(); i++) {
      inbag_counts[oob_sampleIDs[i]] = 0;
    }
  }
}

void BootstrapSampler::bootstrapWithoutReplacementWeighted(std::vector<std::vector<size_t>> &sampleIDs) {

// Use fraction (default 63.21%) of the samples
  size_t num_samples_inbag = (size_t) num_samples * options->get_sample_fraction();
  drawWithoutReplacementWeighted(sampleIDs[0],
                                 num_samples - 1,
                                 num_samples_inbag,
                                 *options->get_case_weights());

// All observation are 0 or 1 times inbag
  inbag_counts.resize(num_samples, 0);
  for (auto &sampleID : sampleIDs[0]) {
    inbag_counts[sampleID] = 1;
  }

  for (size_t s = 0; s < inbag_counts.size(); ++s) {
    if (inbag_counts[s] == 0) {
      oob_sampleIDs.push_back(s);
    }
  }
  num_samples_oob = oob_sampleIDs.size();

  if (!options->get_keep_in_bag()) {
    inbag_counts.clear();
  }
}

void BootstrapSampler::shuffleAndSplit(std::vector<size_t> &first_part,
                                       std::vector<size_t> &second_part,
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

void BootstrapSampler::drawWithoutReplacementSkip(std::vector<size_t> &result,
                                                  size_t max,
                                                  std::vector<size_t> &skip,
                                                  size_t num_samples) {
  if (num_samples < max / 2) {
    drawWithoutReplacementSimple(result, max, skip, num_samples);
  } else {
    drawWithoutReplacementKnuth(result, max, skip, num_samples);
  }
}

void BootstrapSampler::drawWithoutReplacementSimple(std::vector<size_t> &result,
                                                    size_t max,
                                                    std::vector<size_t> &skip,
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
      for (auto &skip_value : skip) {
        if (draw >= skip_value) {
          ++draw;
        }
      }
    } while (temp[draw]);
    temp[draw] = true;
    result.push_back(draw);
  }
}

void BootstrapSampler::drawWithoutReplacementKnuth(std::vector<size_t> &result,
                                                   size_t max,
                                                   std::vector<size_t> &skip,
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
      for (auto &skip_value : skip) {
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

void BootstrapSampler::drawWithoutReplacementWeighted(std::vector<size_t> &result,
                                                      std::vector<size_t> &indices,
                                                      size_t num_samples,
                                                      std::vector<double> &weights) {
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

void BootstrapSampler::drawWithoutReplacementWeighted(std::vector<size_t> &result,
                                                      size_t max_index,
                                                      size_t num_samples,
                                                      std::vector<double> &weights) {
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