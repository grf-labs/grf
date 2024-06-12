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

 #include "relabeling/MultiNoopRelabelingStrategy.h"

 namespace grf {

 MultiNoopRelabelingStrategy::MultiNoopRelabelingStrategy(size_t num_outcomes) :
  num_outcomes(num_outcomes) {}

 bool MultiNoopRelabelingStrategy::relabel(
     const std::vector<size_t>& samples,
     const Data& data,
     Eigen::ArrayXXd& responses_by_sample) const {

   for (size_t sample : samples) {
     responses_by_sample.row(sample) = data.get_outcomes(sample);
   }
   return false;
 }

size_t MultiNoopRelabelingStrategy::get_response_length() const {
  return num_outcomes;
}

 } // namespace grf
