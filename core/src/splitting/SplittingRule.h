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

#ifndef GRF_SPLITTINGRULE_H
#define GRF_SPLITTINGRULE_H

#include <unordered_map>
#include <vector>

class SplittingRule {
public:
  virtual ~SplittingRule() {}
  virtual bool find_best_split(size_t node,
                               const std::vector<size_t>& possible_split_vars,
                               const std::unordered_map<size_t, double>& labels_by_sample,
                               const std::vector<std::vector<size_t>>& samples,
                               std::vector<size_t>& split_vars,
                               std::vector<double>& split_values) = 0;
};

#endif //GRF_SPLITTINGRULE_H


