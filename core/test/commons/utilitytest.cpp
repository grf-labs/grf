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
#include "commons/utility.h"

using namespace grf;

TEST_CASE("Split 0..9 in 1 part", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 0, 9, 1);
  REQUIRE(std::vector<uint>({0, 10}) == test);
}

TEST_CASE("Split 0..7 in 4 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 0, 7, 4);
  REQUIRE(std::vector<uint>({0, 2, 4, 6, 8}) == test);
}

TEST_CASE("Split 2..7 in 2 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 2, 7, 2);
  REQUIRE(std::vector<uint>({2, 5, 8}) == test);
}

TEST_CASE("Split 13..24 in 3 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 13, 24, 3);
  REQUIRE(std::vector<uint>({13, 17, 21, 25}) == test);
}

TEST_CASE("Split 0..6 in 4 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 0, 6, 4);
  REQUIRE(std::vector<uint>({0, 2, 4, 6, 7}) == test);
}

TEST_CASE("Split 2..12 in 5 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 2, 12, 5);
  REQUIRE(std::vector<uint>({2, 5, 7, 9, 11, 13}) == test);
}

TEST_CASE("Split 15..19 in 2 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 15, 19, 2);
  REQUIRE(std::vector<uint>({15, 18, 20}) == test);
}

TEST_CASE("Split 30..35 in 1 parts", "[equalSplit]") {
  std::vector<uint> test;
  split_sequence(test, 30, 35, 1);
  REQUIRE(std::vector<uint>({30, 36}) == test);
}

TEST_CASE("Split 0..2 in 6 parts", "[equalParts]") {
  std::vector<uint> test;
  split_sequence(test, 0, 2, 6);
  REQUIRE(std::vector<uint>({0, 1, 2, 3}) == test);
}

TEST_CASE("Split 0..2 in 4 parts", "[equalParts]") {
  std::vector<uint> test;
  split_sequence(test, 0, 2, 4);
  REQUIRE(std::vector<uint>({0, 1, 2, 3}) == test);
}

TEST_CASE("Split 0..2 in 3 parts", "[equalParts]") {
  std::vector<uint> test;
  split_sequence(test, 0, 2, 3);
  REQUIRE(std::vector<uint>({0, 1, 2, 3}) == test);
}
