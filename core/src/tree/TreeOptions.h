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

namespace grf {

class TreeOptions {
public:
  TreeOptions(uint mtry,
              uint min_node_size,
              bool honesty,
              double honesty_fraction,
              bool prune_empty_leaves,
              double alpha,
              double imbalance_penalty);

  uint get_mtry() const;
  uint get_min_node_size() const;

  /**
  * get_honesty(), get_honesty_fraction(), get_prune_empty_leaves() are all related:
  * if honesty is false, the latter two have no effect, if it is true,
  * get_honesty_fraction is the fraction used to determine splits, and
  * if prune_empty_leaves is true, the resulting honest tree is pruned such that the
  * tree in the estimation sample does not contain any empty leafs. If false,
  * the resulting tree may contain empty leafs, which are skipped at prediction.
  */
  bool get_honesty() const;
  double get_honesty_fraction() const;
  bool get_prune_empty_leaves() const;

  /**
   * The minimum fraction of samples that are allowed to be on either
   * side of each tree split. Splits that are too uneven according to this
   * parameter will not be considered.
   */
  double get_alpha() const;

  /**
   * A tuning parameter that controls how harshly imbalanced splits are penalized.
   */
  double get_imbalance_penalty() const;

private:
  uint mtry;
  uint min_node_size;
  bool honesty;
  double honesty_fraction;
  bool prune_empty_leaves;
  double alpha;
  double imbalance_penalty;
};

} // namespace grf

#endif //GRF_TREEOPTIONS_H
