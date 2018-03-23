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

#include "SamplingOptions.h"

SamplingOptions::SamplingOptions():
    sample_with_replacement(false),
    sample_weights(0) {}

SamplingOptions::SamplingOptions(bool sample_with_replacement):
    sample_with_replacement(sample_with_replacement),
    sample_weights(0) {}

bool SamplingOptions::get_sample_with_replacement() const {
  return sample_with_replacement;
}

const std::vector<double>& SamplingOptions::get_sample_weights() const {
  return sample_weights;
}