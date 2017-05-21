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

#ifndef GRF_PREDICTIONRESULTS_H
#define GRF_PREDICTIONRESULTS_H

#include <cstddef>
#include <vector>

class Prediction {
public:
  Prediction(const std::vector<double>& predictions);
  Prediction(const std::vector<double>& predictions,
             const std::vector<double>& variance_estimates);

  const std::vector<double>& get_predictions() {
    return predictions;
  }

  const std::vector<double>& get_variance_estimates() {
    return variance_estimates;
  }

  bool contains_variance_estimates() {
    return !variance_estimates.empty();
  }

  size_t size() {
    return predictions.size();
  }

private:
  std::vector<double> predictions;
  std::vector<double> variance_estimates;
};


#endif //GRF_PREDICTIONRESULTS_H
