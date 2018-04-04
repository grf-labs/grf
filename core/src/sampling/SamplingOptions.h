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

#ifndef GRF_SAMPLINGOPTIONS_H
#define GRF_SAMPLINGOPTIONS_H


#include <string>
#include <unordered_map>
#include <vector>

#include "commons/globals.h"

class SamplingOptions {
public:
  SamplingOptions();
  SamplingOptions(bool sample_with_replacement, uint samples_per_cluster, std::vector<uint>& clusters);

  bool get_sample_with_replacement() const;
  const std::vector<double>& get_sample_weights() const;

  bool clustering_enabled() const;
  uint get_samples_per_cluster() const;
  std::unordered_map<uint, std::vector<size_t>>& get_cluster_map();
  size_t get_num_clusters() const;

private:
  bool sample_with_replacement;
  std::vector<double> sample_weights;
  uint samples_per_cluster;
  std::unordered_map<uint, std::vector<size_t>> cluster_map;
};

#endif //GRF_SAMPLINGOPTIONS_H
