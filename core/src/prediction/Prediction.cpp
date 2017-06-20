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

Prediction::Prediction(const std::vector<double>& predictions):
  predictions(predictions), variance_estimates(0) {}

Prediction::Prediction(const std::vector<double>& predictions,
                       const std::vector<double>& variance_estimates):
  predictions(predictions), variance_estimates(variance_estimates) {}
