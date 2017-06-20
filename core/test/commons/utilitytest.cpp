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

#include "catch.hpp"
#include "commons/utility.h"

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

TEST_CASE("Read and write a 1D vector", "[fileIO]") {
  std::ofstream outfile;
  outfile.open("testfile1d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<int> expect = std::vector<int>({1, 2, 3, 4, 5, 6});
  write_vector(expect, outfile);
  outfile.close();

  std::ifstream infile;
  infile.open("testfile1d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<int> test;
  read_vector(test, infile);
  infile.close();

  std::remove("testfile1d");

  REQUIRE(expect == test);
}

TEST_CASE("Read and write a 1D vector of doubles", "[fileIO]") {
  std::ofstream outfile;
  outfile.open("testfile1d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<double> expect = std::vector<double>({1.5, 4.5, 10.2, 0.0, -5.9});
  write_vector(expect, outfile);
  outfile.close();

  std::ifstream infile;
  infile.open("testfile1d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<double> test;
  read_vector(test, infile);
  infile.close();

  std::remove("testfile1d");

  REQUIRE(expect == test);
}

TEST_CASE("Read and write a 2D vector of integers", "[fileIO]") {
  std::ofstream outfile;
  outfile.open("testfile2d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<int> expect1 = std::vector<int>({1, 2, 3, 4, 5, 6});
  std::vector<int> expect2 = std::vector<int>({7, 8, 9, 10});
  std::vector<int> expect3 = std::vector<int>({11, 12, 13, 14, 15});
  std::vector<std::vector<int>> expect;
  expect.push_back(expect1);
  expect.push_back(expect2);
  expect.push_back(expect3);

  write_matrix(expect, outfile);
  outfile.close();

  std::ifstream infile;
  infile.open("testfile2d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<std::vector<int>> test;
  read_matrix(test, infile);
  infile.close();

  std::remove("testfile2d");

  REQUIRE(expect == test);
}

TEST_CASE("Read and write a 2D vector of doubles", "[fileIO]") {
  std::ofstream outfile;
  outfile.open("testfile2d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<double> expect1 = std::vector<double>({1.1, 2.4, 3.0, 4.3, 5.9, 6.7});
  std::vector<double> expect2 = std::vector<double>({7.2, 8.1, 9.0, 10.1});
  std::vector<double> expect3 = std::vector<double>({11.3, 12.4, 13.2, 14.7, 15.8});
  std::vector<std::vector<double>> expect;
  expect.push_back(expect1);
  expect.push_back(expect2);
  expect.push_back(expect3);

  write_matrix(expect, outfile);
  outfile.close();

  std::ifstream infile;
  infile.open("testfile2d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<std::vector<double>> test;
  read_matrix<double>(test, infile);
  infile.close();

  std::remove("testfile2d");

  REQUIRE(expect == test);
}

TEST_CASE("Beautify one second", "[beautifyTime]") {
  REQUIRE("1 seconds" == beautify_time(1));
}

TEST_CASE("Beautify seconds", "[beautifyTime]") {
  REQUIRE("30 seconds" == beautify_time(30));
}

TEST_CASE("Beautify one minute", "[beautifyTime]") {
  REQUIRE("1 minute, 0 seconds" == beautify_time(60));
}

TEST_CASE("Beautify minutes", "[beautifyTime]") {
  REQUIRE("38 minutes, 37 seconds" == beautify_time(2317));
}

TEST_CASE("Beautify one hour", "[beautifyTime]") {
  REQUIRE("1 hour, 0 minutes, 0 seconds" == beautify_time(3600));
}

TEST_CASE("Beautify hours", "[beautifyTime]") {
  REQUIRE("3 hours, 44 minutes, 58 seconds" == beautify_time(13498));
}

TEST_CASE("Beautify one day", "[beautifyTime]") {
  REQUIRE("1 day, 0 hours, 0 minutes, 0 seconds" == beautify_time(86400));
}

TEST_CASE("Beautify days", "[beautifyTime]") {
  REQUIRE("3 days, 7 hours, 49 minutes, 5 seconds" == beautify_time(287345));
}

TEST_CASE("Round to next multiple", "[roundToNextMultiple]") {
  REQUIRE(0 == round_to_next_multiple(0, 4));
}

TEST_CASE("Round to next multiple 1", "[roundToNextMultiple]") {
  REQUIRE(4 == round_to_next_multiple(1, 4));
}

TEST_CASE("Round to next multiple 2", "[roundToNextMultiple]") {
  REQUIRE(4 == round_to_next_multiple(2, 4));
}

TEST_CASE("Round to next multiple 3", "[roundToNextMultiple]") {
  REQUIRE(4 == round_to_next_multiple(3, 4));
}

TEST_CASE("Round to next multiple 4", "[roundToNextMultiple]") {
  REQUIRE(4 == round_to_next_multiple(4, 4));
}

TEST_CASE("Round to next multiple 5", "[roundToNextMultiple]") {
  REQUIRE(8 == round_to_next_multiple(5, 4));
}

TEST_CASE("Round to next multiple 6", "[roundToNextMultiple]") {
  REQUIRE(8 == round_to_next_multiple(6, 4));
}

TEST_CASE("Round to next multiple 7", "[roundToNextMultiple]") {
  REQUIRE(8 == round_to_next_multiple(7, 4));
}

TEST_CASE("Round to next multiple 8", "[roundToNextMultiple]") {
  REQUIRE(8 == round_to_next_multiple(8, 4));
}

TEST_CASE("Round to next multiple 9", "[roundToNextMultiple]") {
  REQUIRE(12 == round_to_next_multiple(9, 4));
}

TEST_CASE("Split string 1", "[splitString]") {
  std::string test_string = "abc,def,ghi";
  std::vector<std::string> splitted_string;
  split_string(splitted_string, test_string, ',');
  std::vector<std::string> expect = std::vector<std::string>({"abc", "def", "ghi"});

  REQUIRE(expect == splitted_string);
}

TEST_CASE("Split string 2", "[splitString]") {
  std::string test_string = "abc";
  std::vector<std::string> splitted_string;
  split_string(splitted_string, test_string, ',');
  std::vector<std::string> expect = std::vector<std::string>({"abc"});

  REQUIRE(expect == splitted_string);
}

TEST_CASE("Split string 3", "[splitString]") {
  std::string test_string = "a-b-c";
  std::vector<std::string> splitted_string;
  split_string(splitted_string, test_string, '-');
  std::vector<std::string> expect = std::vector<std::string>({"a", "b", "c"});

  REQUIRE(expect == splitted_string);
}

TEST_CASE("Split string 4", "[splitString]") {
  std::string test_string = "";
  std::vector<std::string> splitted_string;
  split_string(splitted_string, test_string, ',');
  std::vector<std::string> expect = std::vector<std::string>();

  REQUIRE(expect == splitted_string);
}
