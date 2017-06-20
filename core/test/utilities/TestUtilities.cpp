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

#include "TestUtilities.h"

Observations TestUtilities::create_observations(std::vector<double> outcome) {
  std::vector<std::vector<double>> observations_by_type = {outcome};
  return Observations(observations_by_type, outcome.size());
}

Observations TestUtilities::create_observations(std::vector<double> outcome,
                                                std::vector<double> treatment,
                                                std::vector<double> instrument) {
  std::vector<std::vector<double>> observations_by_type = {outcome, treatment, instrument};
  return Observations(observations_by_type, outcome.size());
}
