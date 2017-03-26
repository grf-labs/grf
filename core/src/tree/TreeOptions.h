/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRADIENTFOREST_TREEOPTIONS_H
#define GRADIENTFOREST_TREEOPTIONS_H


#include <string>
#include <vector>

#include "commons/globals.h"

class TreeOptions {
public:
  TreeOptions(uint mtry,
              uint min_node_size,
              const std::vector<double>& split_select_weights,
              const std::vector<size_t>& split_select_varIDs,
              const std::vector<size_t>& deterministic_varIDs,
              const std::vector<size_t>& no_split_variables,
              bool honesty);

  uint get_mtry();
  uint get_min_node_size();
  const std::vector<double>& get_split_select_weights();
  const std::vector<size_t>& get_split_select_varIDs();

  const std::vector<size_t>& get_deterministic_varIDs();
  const std::vector<size_t>& get_no_split_variables();

  bool get_honesty();

private:
  uint mtry;
  uint min_node_size;
  std::vector<double> split_select_weights;
  std::vector<size_t> split_select_varIDs;
  std::vector<size_t> deterministic_varIDs;
  std::vector<size_t> no_split_variables;

  bool honesty;
};

#endif //GRADIENTFOREST_TREEOPTIONS_H
