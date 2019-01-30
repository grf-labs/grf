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

#include "CustomRelabelingStrategy.h"

std::unordered_map<size_t, double> CustomRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data* data) {
  
  std::vector<double> time;
  std::vector<double> status;
  time.reserve(samples.size());
  status.reserve(samples.size());
  for (size_t sample : samples) {
    time.push_back(data->get_outcome(sample));
    status.push_back(data->get_instrument(sample));
  }
  
  std::vector<double> scores = logrankScores(time, status);
  
  std::unordered_map<size_t, double> relabeled_observations;
  for (size_t i = 0; i < samples.size(); ++i) {
    size_t sample = samples[i];
    double score = scores[i];
    relabeled_observations[sample] = score;
  }
  return relabeled_observations;
}

std::vector<double> logrankScores(const std::vector<double>& time, const std::vector<double>& status) {
  size_t n = time.size();
  std::vector<double> scores(n);
  
  // Get order of timepoints
  std::vector<size_t> indices = order(time, false);
  
  // Compute scores
  double cumsum = 0;
  size_t last_unique = -1;
  for (size_t i = 0; i < n; ++i) {
    
    // Continue if next value is the same
    if (i < n - 1 && time[indices[i]] == time[indices[i + 1]]) {
      continue;
    }
    
    // Compute sum and scores for all non-unique values in a row
    for (size_t j = last_unique + 1; j <= i; ++j) {
      cumsum += status[indices[j]] / (n - i);
    }
    for (size_t j = last_unique + 1; j <= i; ++j) {
      scores[indices[j]] = status[indices[j]] - cumsum;
    }
    
    // Save last computed value
    last_unique = i;
  }
  
  return scores;
}
