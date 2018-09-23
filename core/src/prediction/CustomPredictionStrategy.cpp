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
#include "Eigen/Dense"

const std::size_t CustomPredictionStrategy::OUTCOME = 0;

size_t CustomPredictionStrategy::prediction_length() {
  return 1;
}

std::vector<double> CustomPredictionStrategy::predict(size_t sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Observations& observations) {
  return { 0.0 };
}

Eigen::MatrixXd CustomPredictionStrategy::find_M(
        std::unordered_map<size_t, double> weights_by_sampleID,
        size_t sampleID,
        const Observations& observations,
        double lambda){
  Eigen::MatrixXd empty_matrix = Eigen::MatrixXd::Zero (1, 1);
  return empty_matrix;
}

std::vector<double> CustomPredictionStrategy::compute_variance(
        const PredictionValues& leaf_values,
        uint ci_group_size,
        size_t sampleID,
        std::unordered_map<size_t, double> weights_by_sampleID,
        const Observations& observations){
  return { 0.0 };
}

std::vector<double> CustomPredictionStrategy::compute_debiased_error(
        size_t sample,
        const PredictionValues& leaf_values,
        const Observations& observations){
  return { 0.0 };
}

size_t CustomPredictionStrategy::prediction_value_length() {
  return 1;
}