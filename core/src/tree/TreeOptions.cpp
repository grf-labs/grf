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
                         bool honesty,
                         double honesty_fraction,
                         double alpha,
                         double imbalance_penalty):
  mtry(mtry),
  min_node_size(min_node_size),
  honesty(honesty),
  honesty_fraction(honesty_fraction),
  alpha(alpha),
  imbalance_penalty(imbalance_penalty) {}

uint TreeOptions::get_mtry() const  {
  return mtry;
}

uint TreeOptions::get_min_node_size() const  {
  return min_node_size;
}

bool TreeOptions::get_honesty() const {
  return honesty;
}

double TreeOptions::get_honesty_fraction() const {
  return honesty_fraction;
}

double TreeOptions::get_alpha() const {
  return alpha;
}

double TreeOptions::get_imbalance_penalty() const {
  return imbalance_penalty;
}
