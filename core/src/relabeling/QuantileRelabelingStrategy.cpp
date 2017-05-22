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

#include <algorithm>
#include <unordered_map>

#include "relabeling/QuantileRelabelingStrategy.h"

QuantileRelabelingStrategy::QuantileRelabelingStrategy(const std::vector<double>& quantiles) :
    quantiles(quantiles) {}

std::unordered_map<size_t, double> QuantileRelabelingStrategy::relabel(
    const std::vector<size_t>& sampleIDs,
    const Observations& observations) {

  std::vector<double> sorted_outcomes;
  for (auto& sampleID : sampleIDs) {
    sorted_outcomes.push_back(observations.get(Observations::OUTCOME, sampleID));
  }
  std::sort(sorted_outcomes.begin(), sorted_outcomes.end());

  size_t num_samples = sorted_outcomes.size();
  std::vector<double> quantile_cutoffs;

  // Calculate the outcome value cutoffs for each quantile.
  for (auto& quantile : quantiles) {
    size_t outcome_index = (size_t) ceil(num_samples * quantile) - 1;
    quantile_cutoffs.push_back(sorted_outcomes.at(outcome_index));
  }

  // Remove duplicate cutoffs.
  quantile_cutoffs.erase(std::unique(quantile_cutoffs.begin(), quantile_cutoffs.end()),
                         quantile_cutoffs.end());

  // Assign a class to each response based on what quantile it belongs to.
  std::unordered_map<size_t, double> relabeled_observations;
  for (size_t sampleID : sampleIDs) {
    double outcome = observations.get(Observations::OUTCOME, sampleID);
    auto quantile = std::lower_bound(quantile_cutoffs.begin(),
                                     quantile_cutoffs.end(),
                                     outcome);
    long quantile_index = quantile - quantile_cutoffs.begin();
    relabeled_observations[sampleID] = (uint) quantile_index;
  }
  return relabeled_observations;
}
