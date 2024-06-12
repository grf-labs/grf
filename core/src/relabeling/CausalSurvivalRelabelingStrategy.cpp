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

#include "commons/utility.h"
#include "relabeling/CausalSurvivalRelabelingStrategy.h"

namespace grf {

bool CausalSurvivalRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    Eigen::ArrayXXd& responses_by_sample) const {

  // Prepare the relevant averages.
  double numerator_sum = 0;
  double denominator_sum = 0;
  double sum_weight = 0.0;

  for (size_t sample : samples) {
    double sample_weight = data.get_weight(sample);
    numerator_sum += sample_weight * data.get_causal_survival_numerator(sample);
    denominator_sum += sample_weight * data.get_causal_survival_denominator(sample);
    sum_weight += sample_weight;
  }

  if (equal_doubles(denominator_sum, 0.0, 1.0e-10) || std::abs(sum_weight) <= 1e-16) {
    return true;
  }

  double tau = numerator_sum / denominator_sum;

  // Create the new outcomes.
  for (size_t sample : samples) {
    double response = (data.get_causal_survival_numerator(sample) -
      data.get_causal_survival_denominator(sample) * tau) / denominator_sum;
    responses_by_sample(sample, 0) = response;
  }
  return false;
}

} // namespace grf
