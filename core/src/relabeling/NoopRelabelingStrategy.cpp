/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest (grf).

  generalized-random-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  generalized-random-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with generalized-random-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include "relabeling/NoopRelabelingStrategy.h"

std::unordered_map<size_t, double> NoopRelabelingStrategy::relabel(
    const std::vector<size_t>& sampleIDs,
    const Observations& observations) {

  std::unordered_map<size_t, double> relabeled_observations;
  for (size_t sampleID : sampleIDs) {
    double outcome = observations.get(Observations::OUTCOME, sampleID);
    relabeled_observations[sampleID] = outcome;
  }
  return relabeled_observations;
}
