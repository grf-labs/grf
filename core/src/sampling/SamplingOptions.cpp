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

#include "SamplingOptions.h"

#include "commons/globals.h"

SamplingOptions::SamplingOptions():
    sample_with_replacement(false),
    sample_weights(0),
    samples_per_cluster(0),
    cluster_map(0) {}

SamplingOptions::SamplingOptions(bool sample_with_replacement, uint samples_per_cluster, std::vector<uint>& clusters):
        sample_with_replacement(sample_with_replacement),
        sample_weights(0),
        samples_per_cluster(samples_per_cluster) {
  // Create map containing all obs for each cluster. Saves on expense of having to recalculate many times in sampler
  for (size_t s = 0; s < clusters.size(); ++s) {
    size_t obs_cluster = clusters[s];
    this->cluster_map[obs_cluster].push_back(s);
  }
}

bool SamplingOptions::get_sample_with_replacement() const {
  return sample_with_replacement;
}

const std::vector<double>& SamplingOptions::get_sample_weights() const {
  return sample_weights;
}

unsigned int SamplingOptions::get_samples_per_cluster() const {
  return samples_per_cluster;
}

std::unordered_map<uint, std::vector<size_t>>& SamplingOptions::get_cluster_map() {
  return cluster_map;
}

bool SamplingOptions::clustering_enabled() const {
  return !cluster_map.empty();
}

size_t SamplingOptions::get_num_clusters() const {
  return cluster_map.size();
}
