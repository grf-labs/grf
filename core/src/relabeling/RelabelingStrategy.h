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

#ifndef GRF_RELABELINGSTRATEGY_H
#define GRF_RELABELINGSTRATEGY_H


#include "commons/Data.h"
#include "commons/Observations.h"
#include <unordered_map>
#include <vector>

/**
 * Produces a relabelled set of outcomes for a set of training samples. These outcomes
 * will then be used in calculating a standard regression (or classification) split.
 *
 * sampleIDs: the subset of samples to relabel.
 * observations: the list of observations for all training samples.
 *
 * returns: a map from sample ID to a relabelled outcome.
 */
class RelabelingStrategy {
public:
  virtual std::unordered_map<size_t, double> relabel(
      const std::vector<size_t>& sampleIDs,
      const Observations& observations) = 0;
};


#endif //GRF_RELABELINGSTRATEGY_H
