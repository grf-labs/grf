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

#include "relabeling/NoopRelabelingStrategy.h"

std::unordered_map<size_t, double> NoopRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data* data) {

  std::unordered_map<size_t, double> relabeled_observations;
  for (size_t sample : samples) {
    double outcome = data->get_outcome(sample);
    relabeled_observations[sample] = outcome;
  }
  return relabeled_observations;
}
