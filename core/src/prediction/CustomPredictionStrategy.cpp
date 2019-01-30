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

#include "CustomPredictionStrategy.h"
 
 CustomPredictionStrategy::CustomPredictionStrategy(std::vector<double> timepoints):
   timepoints(timepoints) {
 };

size_t CustomPredictionStrategy::prediction_length() {
  return timepoints.size();
}

std::vector<double> CustomPredictionStrategy::predict(size_t sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Data* train_data,
    const Data* data) {
  
  size_t n = weights_by_sample.size();
  
  std::vector<double> time;
  std::vector<double> status;
  std::vector<double> weights;
  time.reserve(n);
  status.reserve(n);
  weights.reserve(n);
  for (auto it = weights_by_sample.begin(); it != weights_by_sample.end(); it++) {
    size_t sample = it->first;
    
    time.push_back(train_data->get_outcome(sample));
    status.push_back(train_data->get_instrument(sample));
    weights.push_back(it->second * (double) n);
  }
  
  return weighted_kaplan_meier(time, status, weights);
}

std::vector<double> CustomPredictionStrategy::compute_variance(
    size_t sample,
    std::vector<std::vector<size_t>> samples_by_tree,
    std::unordered_map<size_t, double> weights_by_sampleID,
    const Data* train_data,
    const Data* data,
    size_t ci_group_size){
  return { 0.0 };
}

std::vector<double> CustomPredictionStrategy::weighted_kaplan_meier(std::vector<double> time, std::vector<double> status,
                                          std::vector<double> weights) {
  
  size_t n = time.size();             
  size_t num_timepoints = timepoints.size();            
  
  // Compute at risk and events for each timepoint
  std::vector<double> at_risk(num_timepoints);
  std::vector<double> death(num_timepoints);
  for (size_t i = 0; i < n; ++i) {
      double survival_time = time[i];
      for(size_t j = 0; j < num_timepoints; ++j) {
        if (survival_time >= timepoints[j]) {
          at_risk[j] += weights[i];
        } else if (status[i] == 0) {
          break;
        } else if (status[i] == 1) {
          at_risk[j] += weights[i];
        }
        
        if (status[i] == 1 && survival_time <= timepoints[j]) {
          death[j] += weights[i];
          break;
        }
      }
  } 
  
  // Compute KM estimator 
  std::vector<double> survival(num_timepoints);
  if (at_risk[0] > 0) {
    survival[0] = 1 - death[0]/at_risk[0];
  }
  for (size_t j = 1; j < num_timepoints; ++j) {
    survival[j] = survival[j-1];
    if (at_risk[j] > 0) {
      survival[j] *= (1 - death[j]/at_risk[j]);
    }
  }
  
  return survival;
}
