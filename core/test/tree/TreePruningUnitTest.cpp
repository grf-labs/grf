/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include "catch.hpp"
#include "tree/Tree.h"

using namespace grf;

TEST_CASE("pruning behaves as expected", "[tree, unit]") {
  /*
   * This test case starts with the following tree structure, with all leaves
   * empty except nodes 6 and 7:
   *
   *             0
   *           /   \
   *          1     2
   *        /   \
   *       3     4
   *      / \   / \
   *     5   6 7   8
   *              / \
   *             9  10
   *
   * Pruning should produce this tree (note the root node ID has changed):
   *
   *          1
   *         / \
   *        6   7
   */

  std::vector<std::vector<size_t>> child_nodes =
      {{1, 3, 0, 5, 7, 0, 0, 0, 9, 0, 0}, {2, 4, 0, 6, 8, 0, 0, 0, 10, 0, 0}};
  std::vector<std::vector<size_t>> leaf_nodes = {
      {{}, {}, {}, {}, {}, {}, {42, 43}, {44}, {}, {}, {}}};
  Tree tree(0, child_nodes, leaf_nodes, {0}, {0}, {0}, {true}, PredictionValues());

  tree.honesty_prune_leaves();

  std::vector<std::vector<size_t>> expected_child_nodes =
      {{0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  const std::vector<std::vector<size_t>>& actual_child_nodes = tree.get_child_nodes();

  REQUIRE(tree.get_root_node() == 1);
  for (size_t node = 0; node < leaf_nodes.size(); node++) {
    REQUIRE(actual_child_nodes[0][node] == expected_child_nodes[0][node]);
    REQUIRE(actual_child_nodes[1][node] == expected_child_nodes[1][node]);
  }
}

TEST_CASE("pruning is idempotent", "[tree, unit]") {
  std::vector<std::vector<size_t>> child_nodes =
      {{1, 3, 0, 5, 7, 0, 0, 0, 9, 0, 0}, {2, 4, 0, 6, 8, 0, 0, 0, 10, 0, 0}};
  std::vector<std::vector<size_t>> leaf_nodes = {
      {{}, {}, {}, {}, {}, {}, {42, 43}, {44}, {}, {}, {}}};
  Tree tree(0, child_nodes, leaf_nodes, {0}, {0}, {0}, {true}, PredictionValues());

  tree.honesty_prune_leaves();
  const std::vector<std::vector<size_t>>& expected_child_nodes = tree.get_child_nodes();

  tree.honesty_prune_leaves();
  const std::vector<std::vector<size_t>>& actual_child_nodes = tree.get_child_nodes();

  for (size_t node = 0; node < leaf_nodes.size(); node++) {
    REQUIRE(actual_child_nodes[0][node] == expected_child_nodes[0][node]);
    REQUIRE(actual_child_nodes[1][node] == expected_child_nodes[1][node]);
  }
}
