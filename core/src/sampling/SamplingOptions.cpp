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
#include <unordered_map>
#include "commons/globals.h"

namespace grf {

SamplingOptions::SamplingOptions():
    sample_weights(0),
    num_samples_per_cluster(0),
    clusters(0) {}

SamplingOptions::SamplingOptions(uint samples_per_cluster,
                                 const std::vector<size_t>& sample_clusters):
    sample_weights(0),
    num_samples_per_cluster(samples_per_cluster) {

  // Map the provided clusters to IDs in the range 0 ... num_clusters.
  std::unordered_map<size_t, size_t> cluster_ids;
  for (size_t cluster : sample_clusters) {
    if (cluster_ids.find(cluster) == cluster_ids.end()) {
      size_t cluster_id = cluster_ids.size();
      cluster_ids[cluster] = cluster_id;
    }
  }

  // Populate the index of each cluster ID with the samples it contains.
  clusters = std::vector<std::vector<size_t>>(cluster_ids.size());
  for (size_t sample = 0; sample < sample_clusters.size(); sample++) {
    size_t cluster = sample_clusters.at(sample);
    size_t cluster_id = cluster_ids.at(cluster);
    clusters[cluster_id].push_back(sample);
  }
}

const std::vector<double>& SamplingOptions::get_sample_weights() const {
  return sample_weights;
}


unsigned int SamplingOptions::get_samples_per_cluster() const {
  return num_samples_per_cluster;
}

const std::vector<std::vector<size_t>>& SamplingOptions::get_clusters() const {
  return clusters;
}

} // namespace grf
