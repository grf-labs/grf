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

#include <map>
#include <unordered_set>
#include <fstream>
#include <random>

#include "catch.hpp"
#include "serialization/ForestSerializer.h"

TEST_CASE("trees serialize and deserialize correctly", "[treeSerialization]") {
  std::shared_ptr<Tree> original_tree(new Tree(0,
      {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
      {{0, 1}, {5, 1}, {4, 6, 9}},
      {10, 20, 30, 40},
      {0.5, 0.75, 0.9, 1.1, 1.2},
      {3, 4, 5, 9, 10, 11},
      PredictionValues({{0, 0, 1}}, 3)));

  TreeSerializer tree_serializer;
  std::stringstream stream;
  tree_serializer.serialize(stream, original_tree);
  std::shared_ptr<Tree> tree = tree_serializer.deserialize(stream);

  REQUIRE(tree->get_leaf_samples().size() == original_tree->get_leaf_samples().size());
  REQUIRE(tree->get_child_nodes().size() == original_tree->get_child_nodes().size());
  REQUIRE(tree->get_split_vars().size() == original_tree->get_split_vars().size());
  REQUIRE(tree->get_split_values().size() == original_tree->get_split_values().size());
  REQUIRE(tree->get_drawn_samples().size() == original_tree->get_drawn_samples().size());
}

TEST_CASE("forests serialize and deserialize correctly", "[forestSerialization]") {
  std::vector<std::shared_ptr<Tree>> trees;
  trees.push_back(std::shared_ptr<Tree>(new Tree(0, {{0}}, {{0}}, {0}, {0}, {0}, PredictionValues())));
  trees.push_back(std::shared_ptr<Tree>(new Tree(0, {{1}}, {{1}}, {1}, {1}, {1}, PredictionValues())));

  ForestSerializer forest_serializer;
  std::stringstream stream;
  Forest original_forest(trees, 3, 5);

  forest_serializer.serialize(stream, original_forest);
  Forest forest = forest_serializer.deserialize(stream);

  REQUIRE(forest.get_trees().size() == original_forest.get_trees().size());
  REQUIRE(forest.get_num_variables() == original_forest.get_num_variables());
  REQUIRE(forest.get_ci_group_size() == original_forest.get_ci_group_size());
}
