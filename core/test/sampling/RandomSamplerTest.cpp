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
#include <fstream>
#include <map>
#include <random>
#include <unordered_set>

#include "catch.hpp"
#include "commons/DefaultData.h"
#include "sampling/RandomSampler.h"

size_t absolute_difference(size_t first, size_t second) {
  return first >= second ? first - second : second - first;
}

TEST_CASE("Draw without replacement 1", "[drawWithoutReplacement]") {
  std::vector<size_t> result;
  std::random_device random_device;
  std::map<size_t, uint> counts;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);

  size_t max = 9;
  std::set<size_t> skip = {7};
  size_t num_samples = 4;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    sampler.draw_without_replacement_skip(result, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    size_t difference = absolute_difference(expected_count, c.second);
    REQUIRE(difference < expected_count * 0.05);
  }
  REQUIRE(0 == counts[*skip.begin()]);
}

TEST_CASE("Draw without replacement 2", "[drawWithoutReplacement]") {
  std::vector<size_t> result;
  std::random_device random_device;
  std::map<size_t, uint> counts;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);

  size_t max = 9;
  std::set<size_t> skip = {0};
  size_t num_samples = 4;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    sampler.draw_without_replacement_skip(result, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    size_t difference = absolute_difference(expected_count, c.second);
    REQUIRE(difference < expected_count * 0.05);
  }
  REQUIRE(0 == counts[*skip.begin()]);
}

TEST_CASE("Draw without replacement 3", "[drawWithoutReplacement]") {
  std::vector<size_t> result;
  std::random_device random_device;
  std::map<size_t, uint> counts;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options); 

  size_t max = 9;
  std::set<size_t> skip = {9};
  size_t num_samples = 4;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    sampler.draw_without_replacement_skip(result, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    size_t difference = absolute_difference(expected_count, c.second);
    REQUIRE(difference < expected_count * 0.05);
  }
  REQUIRE(0 == counts[*skip.begin()]);
}

TEST_CASE("Draw without replacement 4", "[drawWithoutReplacement]") {
  std::vector<size_t> result;
  std::random_device random_device;
  std::map<size_t, uint> counts;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);
  
  size_t max = 1000;
  std::set<size_t> skip = {7};
  size_t num_samples = 50;
  size_t num_replicates = 100000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    sampler.draw_without_replacement_skip(result, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 10%
  for (auto& c : counts) {
    size_t difference = absolute_difference(expected_count, c.second);
    REQUIRE(difference < expected_count * 0.10);
  }
  REQUIRE(0 == counts[*skip.begin()]);
}

TEST_CASE("Draw without replacement 5", "[drawWithoutReplacement]") {
  std::vector<size_t> result;
  std::random_device random_device;
  std::map<size_t, uint> counts;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);

  size_t max = 1000;
  std::set<size_t> skip = {7};
  size_t num_samples = 500;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    sampler.draw_without_replacement_skip(result, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    size_t difference = absolute_difference(expected_count, c.second);
    REQUIRE(difference < expected_count * 0.05);
  }
  REQUIRE(0 == counts[*skip.begin()]);
}


TEST_CASE("Shuffle and split 1", "[shuffleAndSplit]") {
  std::random_device random_device;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);

  std::vector<size_t> samples;

  sampler.shuffle_and_split(samples, 10, 3);

  REQUIRE(3 == samples.size());
}

TEST_CASE("Shuffle and split 2", "[shuffleAndSplit]") {
  std::random_device random_device;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);

  std::vector<size_t> samples;

  sampler.shuffle_and_split(samples, 100, 63);

  REQUIRE(63 == samples.size());
}

TEST_CASE("Shuffle and split 3", "[shuffleAndSplit]") {
  std::random_device random_device;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);

  std::vector<size_t> samples;

  sampler.shuffle_and_split(samples, 1, 1);

  REQUIRE(1 == samples.size());
}

TEST_CASE("Shuffle and split 4", "[shuffleAndSplit]") {
  std::random_device random_device;

  std::vector<uint> dummy_clusters;
  SamplingOptions sampling_options(true, 0, dummy_clusters);
  RandomSampler sampler(random_device(), sampling_options);
  
  std::vector<size_t> samples;

  sampler.shuffle_and_split(samples, 3, 0);

  REQUIRE(0 == samples.size());
}

TEST_CASE("sample multilevel 1", "[sampleMultilevel]") {
  std::random_device random_device;
  DefaultData data(NULL, 0, 0);
  uint samples_per_cluster = 3;
  std::vector<uint> clusters = {0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 2, 2, 2, 2, 0, 3, 3, 3, 2, 3};
  size_t num_clusters = 4;
  SamplingOptions sampling_options(true, samples_per_cluster, clusters);
  RandomSampler sampler(random_device(), sampling_options);

  std::vector<size_t> sampled_clusters;
  sampler.sample_for_ci(&data, .5, sampled_clusters);

  size_t expected_num_sampled_clusters = (size_t) std::ceil(.5 * num_clusters);
  REQUIRE(expected_num_sampled_clusters == sampled_clusters.size());

  std::vector<size_t> subsampled_clusters;
  std::vector<size_t> oob_subsampled_clusters;
  sampler.subsample(sampled_clusters, 0.5, subsampled_clusters, oob_subsampled_clusters);

  std::vector<size_t> samples;
  std::vector<size_t> oob_sample;
  sampler.sample_from_clusters(subsampled_clusters, samples);
  sampler.get_oob_from_clusters(oob_subsampled_clusters, oob_sample);

  std::unordered_set<size_t> samples_set(samples.begin(), samples.end());
  std::unordered_set<size_t> oob_sample_set(oob_sample.begin(), oob_sample.end());
  size_t expected_num_subsampled_clusters = (size_t) std::ceil(.5 * expected_num_sampled_clusters);
  size_t expected_sample_size = expected_num_subsampled_clusters * samples_per_cluster;

  // Check that sample is of correct size
  REQUIRE(expected_sample_size == samples.size());

  // Check that sample and oob_sample don't intersect
  std::vector<size_t> intersection;
  std::set_intersection(samples_set.begin(), samples_set.end(),
                        oob_sample_set.begin(), oob_sample_set.end(),
                        std::back_inserter(intersection));
  REQUIRE(intersection.empty());

  // Check that sample is from appropriate clusters
  std::set<size_t> expected_subsampled_clusters(subsampled_clusters.begin(), subsampled_clusters.end());
  std::set<size_t> actual_subsampled_clusters;
  for (auto const& i: samples) {
    size_t cluster = clusters[i];
      actual_subsampled_clusters.insert(cluster);
  }
  REQUIRE(actual_subsampled_clusters == expected_subsampled_clusters);

  // Check that oob sample is from appropriate clusters
  std::set<size_t> expected_oob_subsampled_clusters(oob_subsampled_clusters.begin(), oob_subsampled_clusters.end());
  std::set<size_t> actual_oob_subsampled_clusters;
  for (auto const& i: oob_sample) {
      size_t cluster = clusters[i];
      actual_oob_subsampled_clusters.insert(cluster);
  }
  REQUIRE(actual_oob_subsampled_clusters == expected_oob_subsampled_clusters);
}
