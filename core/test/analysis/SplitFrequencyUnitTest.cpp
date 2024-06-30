/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

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

#include "catch.hpp"

#include "analysis/SplitFrequencyComputer.h"
#include "forest/Forest.h"

using namespace grf;

TEST_CASE("split frequency computation works as expected", "[analysis, unit]") {
  /*
   * This test case creates the following two trees. The values below represent
   * the variable on which each node was split.
   *
   *             1
   *           /   \
   *          3     x
   *        /   \
   *       4     x
   *      / \
   *     x   x
   *
   *             0
   *           /   \
   *          4     3
   *         / \   / \
   *        x   x x   x
   */
  std::vector<size_t> first_split_vars = {1, 3, 0, 4, 0, 0, 0};
  std::vector<std::vector<size_t>> first_child_nodes =
      {{1, 3, 0, 5, 0, 0, 0}, {2, 4, 0, 6, 0, 0, 0}};

  std::vector<size_t> second_split_vars = {0, 4, 3, 0, 0, 0, 0};
  std::vector<std::vector<size_t>> second_child_nodes =
      {{1, 3, 5, 0, 0, 0, 0}, {2, 4, 6, 0, 0, 0, 0}};

  std::vector<std::vector<size_t>> expected_variable_frequencies = {
      {1, 1, 0, 0, 0}, // depth 1
      {0, 0, 0, 2, 1}, // depth 2
      {0, 0, 0, 0, 1}}; // depth 3

  std::vector<std::unique_ptr<Tree>> trees;
  trees.emplace_back(new Tree(0, first_child_nodes, {{0}}, first_split_vars, {0}, {0}, {true}, PredictionValues()));
  trees.emplace_back(new Tree(0, second_child_nodes, {{1}}, second_split_vars, {1}, {1}, {true}, PredictionValues()));

  size_t num_variables = 5;
  size_t ci_group_size = 2;
  Forest forest(trees, num_variables, ci_group_size);

  SplitFrequencyComputer computer;
  size_t max_depth = 3;
  std::vector<std::vector<size_t>> actual_variable_frequencies = computer.compute(forest, max_depth);

  REQUIRE (actual_variable_frequencies.size() <= max_depth);
  REQUIRE(actual_variable_frequencies.size() == expected_variable_frequencies.size());
  for (size_t var = 0; var < actual_variable_frequencies.size(); var++) {
    const std::vector<size_t>& actual_frequencies = actual_variable_frequencies.at(var);
    const std::vector<size_t>& expected_frequencies = expected_variable_frequencies.at(var);

    REQUIRE(actual_frequencies.size() == num_variables);
    REQUIRE(actual_frequencies.size() == expected_frequencies.size());
    for (size_t depth = 0; depth < actual_frequencies.size(); depth++) {
      REQUIRE(actual_frequencies[depth] == expected_frequencies[depth]);
    }
  }
}

TEST_CASE("split frequency computation respects max depth", "[analysis, unit]") {
  std::vector<size_t> split_vars = {1, 3, 0, 4, 0, 0, 0};
  std::vector<std::vector<size_t>> child_nodes =
      {{1, 3, 0, 5, 0, 0, 0}, {2, 4, 0, 6, 0, 0, 0}};

  std::vector<std::vector<size_t>> expected_variable_frequencies = {
      {1, 1, 0, 0, 0}, // depth 1
      {0, 0, 0, 2, 1}}; // depth 2

  std::vector<std::unique_ptr<Tree>> trees;
  trees.emplace_back(new Tree(0, child_nodes, {{0}}, split_vars, {0}, {0}, {true}, PredictionValues()));

  size_t num_variables = 5;
  size_t ci_group_size = 2;
  Forest forest(trees, num_variables, ci_group_size);

  SplitFrequencyComputer computer;
  size_t max_depth = 2;
  std::vector<std::vector<size_t>> actual_variable_frequencies = computer.compute(forest, max_depth);

  REQUIRE (actual_variable_frequencies.size() <= max_depth);
}
