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

#include "prediction/Prediction.h"

namespace grf {

Prediction::Prediction(const std::vector<double>& predictions):
  predictions(predictions),
  variance_estimates(0),
  error_estimates(0),
  excess_error_estimates(0) {}

Prediction::Prediction(const std::vector<double>& predictions,
                       const std::vector<double>& variance_estimates,
                       const std::vector<double>& error_estimates,
                       const std::vector<double>& excess_error_estimates):
  predictions(predictions),
  variance_estimates(variance_estimates),
  error_estimates(error_estimates),
  excess_error_estimates(excess_error_estimates) {}

const std::vector<double>& Prediction::get_predictions() const {
  return predictions;
}

const std::vector<double>& Prediction::get_variance_estimates() const {
  return variance_estimates;
}

const std::vector<double>& Prediction::get_error_estimates() const {
  return error_estimates;
}

const std::vector<double>& Prediction::get_excess_error_estimates() const {
  return excess_error_estimates;
}

const bool Prediction::contains_variance_estimates() const {
  return !variance_estimates.empty();
}

const bool Prediction::contains_error_estimates() const {
  return !error_estimates.empty();
}

const size_t Prediction::size() const {
  return predictions.size();
}

} // namespace grf
