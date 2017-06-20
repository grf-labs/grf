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

#ifndef GRF_TREEOPTIONS_H
#define GRF_TREEOPTIONS_H


#include <set>
#include <string>
#include <vector>

#include "commons/globals.h"

class TreeOptions {
public:
  TreeOptions(uint mtry,
              uint min_node_size,
              const std::vector<double>& split_select_weights,
              const std::vector<size_t>& split_select_vars,
              const std::vector<size_t>& deterministic_vars,
              const std::set<size_t>& no_split_variables,
              bool honesty);

  uint get_mtry();
  uint get_min_node_size();
  const std::vector<double>& get_split_select_weights();
  const std::vector<size_t>& get_split_select_vars();

  const std::vector<size_t>& get_deterministic_vars();
  const std::set<size_t>& get_no_split_variables();

  bool get_honesty();

private:
  uint mtry;
  uint min_node_size;
  std::vector<double> split_select_weights;
  std::vector<size_t> split_select_vars;
  std::vector<size_t> deterministic_vars;
  std::set<size_t> no_split_variables;

  bool honesty;
};

#endif //GRF_TREEOPTIONS_H
