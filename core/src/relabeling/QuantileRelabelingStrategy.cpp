/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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
#include "relabeling/QuantileRelabelingStrategy.h"

namespace grf {

QuantileRelabelingStrategy::QuantileRelabelingStrategy(const std::vector<double>& quantiles) :
    quantiles(quantiles) {}

bool QuantileRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    Eigen::ArrayXXd& responses_by_sample) const {

  std::vector<double> sorted_outcomes(samples.size());
  for (size_t i = 0; i < samples.size(); i++) {
    size_t sample = samples[i];
    sorted_outcomes[i] = data.get_outcome(sample);
  }
  std::sort(sorted_outcomes.begin(), sorted_outcomes.end());

  size_t num_samples = sorted_outcomes.size();
  std::vector<double> quantile_cutoffs;

  // Calculate the outcome value cutoffs for each quantile.
  for (double quantile : quantiles) {
    size_t outcome_index = (size_t) std::ceil(num_samples * quantile) - 1;
    quantile_cutoffs.push_back(sorted_outcomes[outcome_index]);
  }

  // Remove duplicate cutoffs.
  quantile_cutoffs.erase(std::unique(quantile_cutoffs.begin(), quantile_cutoffs.end()),
                         quantile_cutoffs.end());

  // Assign a class to each response based on what quantile it belongs to.
  for (size_t sample : samples) {
    double outcome = data.get_outcome(sample);
    auto quantile = std::lower_bound(quantile_cutoffs.begin(),
                                     quantile_cutoffs.end(),
                                     outcome);
    long quantile_index = static_cast<long>(quantile - quantile_cutoffs.begin());
    responses_by_sample(sample, 0) = static_cast<uint>(quantile_index);
  }
  return false;
}

} // namespace grf
