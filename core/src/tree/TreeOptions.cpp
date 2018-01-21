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

#include "tree/TreeOptions.h"


TreeOptions::TreeOptions(uint mtry,
                         uint min_node_size,
                         bool honesty):
  mtry(mtry),
  min_node_size(min_node_size),
  honesty(honesty) {}

uint TreeOptions::get_mtry() const  {
  return mtry;
}

uint TreeOptions::get_min_node_size() const  {
  return min_node_size;
}

bool TreeOptions::get_honesty() const {
  return honesty;
}

void TreeOptions::set_min_node_size(uint min_node_size) {
  this->min_node_size = min_node_size;
}
