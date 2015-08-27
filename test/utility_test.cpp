#include <map>
#include <unordered_set>
#include <fstream>

#include "gtest/gtest.h"
#include "utility.h"

// Split 0..9 in 1 part
TEST(equalSplit, onePart) {
  std::vector<uint> test;
  equalSplit(test, 0, 9, 1);
  EXPECT_EQ(std::vector<uint>( { 0, 10 }), test);
}

// Split 0..7 in 4 parts
TEST(equalSplit, perfectSplit0) {
  std::vector<uint> test;
  equalSplit(test, 0, 7, 4);
  EXPECT_EQ(std::vector<uint>( { 0, 2, 4, 6, 8 }), test);
}

// Split 2..7 in 2 parts
TEST(equalSplit, perfectSplit2) {
  std::vector<uint> test;
  equalSplit(test, 2, 7, 2);
  EXPECT_EQ(std::vector<uint>( { 2, 5, 8 }), test);
}

// Split 13..24 in 3 parts
TEST(equalSplit, perfectSplit13) {
  std::vector<uint> test;
  equalSplit(test, 13, 24, 3);
  EXPECT_EQ(std::vector<uint>( { 13, 17, 21, 25 }), test);
}

// Split 0..6 in 4 parts
TEST(equalSplit, nonPerfectSplit0) {
  std::vector<uint> test;
  equalSplit(test, 0, 6, 4);
  EXPECT_EQ(std::vector<uint>( { 0, 2, 4, 6, 7 }), test);
}

// Split 2..12 in 5 parts
TEST(equalSplit, nonPerfectSplit2) {
  std::vector<uint> test;
  equalSplit(test, 2, 12, 5);
  EXPECT_EQ(std::vector<uint>( { 2, 5, 7, 9, 11, 13 }), test);
}

// Split 15..19 in 2 parts
TEST(equalSplit, nonPerfectSplit15) {
  std::vector<uint> test;
  equalSplit(test, 15, 19, 2);
  EXPECT_EQ(std::vector<uint>( { 15, 18, 20 }), test);
}

// Split 30..35 in 1 parts
TEST(equalSplit, nonPerfectSplit30) {
  std::vector<uint> test;
  equalSplit(test, 30, 35, 1);
  EXPECT_EQ(std::vector<uint>( { 30, 36 }), test);
}

// Split 0..2 in 6 parts
// Result should contain only 3 parts
TEST(equalSplit, moreParts1) {
  std::vector<uint> test;
  equalSplit(test, 0, 2, 6);
  EXPECT_EQ(std::vector<uint>( { 0, 1, 2, 3 }), test);
}

// Split 0..2 in 4 parts
// Result should contain only 3 parts
TEST(equalSplit, moreParts2) {
  std::vector<uint> test;
  equalSplit(test, 0, 2, 4);
  EXPECT_EQ(std::vector<uint>( { 0, 1, 2, 3 }), test);
}

// Split 0..2 in 3 parts
// Result should contain only 3 parts
TEST(equalSplit, moreParts3) {
  std::vector<uint> test;
  equalSplit(test, 0, 2, 3);
  EXPECT_EQ(std::vector<uint>( { 0, 1, 2, 3 }), test);
}

TEST(readWrite1D, int1) {
  std::ofstream outfile = std::ofstream();
  outfile.open("testfile1d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<int> expect = std::vector<int>( { 1, 2, 3, 4, 5, 6 });
  saveVector1D(expect, outfile);
  outfile.close();

  std::ifstream infile = std::ifstream();
  infile.open("testfile1d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<int> test;
  readVector1D(test, infile);
  infile.close();

  EXPECT_EQ(expect, test);
}

TEST(readWrite1D, double1) {
  std::ofstream outfile = std::ofstream();
  outfile.open("testfile1d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<double> expect = std::vector<double>( { 1.5, 4.5, 10.2, 0.0, -5.9 });
  saveVector1D(expect, outfile);
  outfile.close();

  std::ifstream infile = std::ifstream();
  infile.open("testfile1d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<double> test;
  readVector1D(test, infile);
  infile.close();

  EXPECT_EQ(expect, test);
}

TEST(readWrite2D, int1) {
  std::ofstream outfile = std::ofstream();
  outfile.open("testfile2d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<int> expect1 = std::vector<int>( { 1, 2, 3, 4, 5, 6 });
  std::vector<int> expect2 = std::vector<int>( { 7, 8, 9, 10 });
  std::vector<int> expect3 = std::vector<int>( { 11, 12, 13, 14, 15 });
  std::vector<std::vector<int>> expect;
  expect.push_back(expect1);
  expect.push_back(expect2);
  expect.push_back(expect3);

  saveVector2D(expect, outfile);
  outfile.close();

  std::ifstream infile = std::ifstream();
  infile.open("testfile2d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<std::vector<int>> test;
  readVector2D(test, infile);
  infile.close();

  EXPECT_EQ(expect, test);
}

TEST(readWrite2D, double1) {
  std::ofstream outfile = std::ofstream();
  outfile.open("testfile2d", std::ios::binary);
  if (!outfile.good()) {
    throw std::runtime_error("Could not write to output file: ");
  }
  std::vector<double> expect1 = std::vector<double>( { 1.1, 2.4, 3.0, 4.3, 5.9, 6.7 });
  std::vector<double> expect2 = std::vector<double>( { 7.2, 8.1, 9.0, 10.1 });
  std::vector<double> expect3 = std::vector<double>( { 11.3, 12.4, 13.2, 14.7, 15.8 });
  std::vector<std::vector<double>> expect;
  expect.push_back(expect1);
  expect.push_back(expect2);
  expect.push_back(expect3);

  saveVector2D(expect, outfile);
  outfile.close();

  std::ifstream infile = std::ifstream();
  infile.open("testfile2d", std::ios::binary);
  if (!infile.good()) {
    throw std::runtime_error("Could not read from input file: ");
  }
  std::vector<std::vector<double>> test;
  readVector2D<double>(test, infile);
  infile.close();

  EXPECT_EQ(expect, test);
}

TEST(drawWithoutReplacementSkip, small_small1) {

  std::vector<size_t> result;
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());
  std::map<size_t, uint> counts;

  size_t max = 9;
  std::vector<size_t> skip = std::vector<size_t>( { 7 });
  size_t num_samples = 4;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    drawWithoutReplacementSkip(result, random_number_generator, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    EXPECT_NEAR(expected_count, c.second, expected_count * 0.05);
  }
  EXPECT_EQ(0, counts[skip[0]]);
}

TEST(drawWithoutReplacementSkip, small_small2) {

  std::vector<size_t> result;
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());
  std::map<size_t, uint> counts;

  size_t max = 9;
  std::vector<size_t> skip = std::vector<size_t>( { 0 });
  size_t num_samples = 4;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    drawWithoutReplacementSkip(result, random_number_generator, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    EXPECT_NEAR(expected_count, c.second, expected_count * 0.05);
  }
  EXPECT_EQ(0, counts[skip[0]]);
}

TEST(drawWithoutReplacementSkip, small_small3) {

  std::vector<size_t> result;
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());
  std::map<size_t, uint> counts;

  size_t max = 9;
  std::vector<size_t> skip = std::vector<size_t>( { 9 });
  size_t num_samples = 4;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    drawWithoutReplacementSkip(result, random_number_generator, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    EXPECT_NEAR(expected_count, c.second, expected_count * 0.05);
  }
  EXPECT_EQ(0, counts[skip[0]]);
}

TEST(drawWithoutReplacementSkip, small_large1) {

  std::vector<size_t> result;
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());
  std::map<size_t, uint> counts;

  size_t max = 1000;
  std::vector<size_t> skip = std::vector<size_t>( { 7 });
  size_t num_samples = 50;
  size_t num_replicates = 100000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    drawWithoutReplacementSkip(result, random_number_generator, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 10%
  for (auto& c : counts) {
    EXPECT_NEAR(expected_count, c.second, expected_count * 0.1);
  }
  EXPECT_EQ(0, counts[skip[0]]);
}

TEST(drawWithoutReplacementSkip, large_large1) {

  std::vector<size_t> result;
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());
  std::map<size_t, uint> counts;

  size_t max = 1000;
  std::vector<size_t> skip = std::vector<size_t>( { 7 });
  size_t num_samples = 500;
  size_t num_replicates = 10000;

  size_t expected_count = num_samples * num_replicates / max;

  for (size_t i = 0; i < num_replicates; ++i) {
    result.clear();
    drawWithoutReplacementSkip(result, random_number_generator, max + 1, skip, num_samples);
    for (auto& idx : result) {
      ++counts[idx];
    }
  }

  // Check if counts are expected +- 5%
  for (auto& c : counts) {
    EXPECT_NEAR(expected_count, c.second, expected_count * 0.05);
  }
  EXPECT_EQ(0, counts[skip[0]]);
}

TEST(mostFrequentClass, notEqual1) {
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<uint> class_count = std::vector<uint>( { 0, 4, 7, 3, 2, 1, 8 });

  EXPECT_EQ(6, mostFrequentClass(class_count, random_number_generator));
}

TEST(mostFrequentClass, notEqual2) {
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<uint> class_count = std::vector<uint>( { 5, 4, 3, 2, 1 });

  EXPECT_EQ(0, mostFrequentClass(class_count, random_number_generator));
}

TEST(mostFrequentClass, equal1) {
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<uint> class_count = std::vector<uint>( { 5, 5, 5, 5 });

  EXPECT_LE(0, mostFrequentClass(class_count, random_number_generator));
  EXPECT_GE(3, mostFrequentClass(class_count, random_number_generator));
}

TEST(mostFrequentClass, equal2) {
  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<uint> class_count = std::vector<uint>( { 4, 5, 5, 4 });

  EXPECT_LE(1, mostFrequentClass(class_count, random_number_generator));
  EXPECT_GE(2, mostFrequentClass(class_count, random_number_generator));
}

TEST(mostFrequentValue, notEqual1) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::unordered_map<double, size_t> class_count;
  class_count[1] = 5;
  class_count[2] = 7;
  class_count[3] = 10;

  EXPECT_EQ(3, mostFrequentValue(class_count, random_number_generator));
}

TEST(mostFrequentValue, notEqual2) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::unordered_map<double, size_t> class_count;
  class_count[10.1] = 15;
  class_count[2.5] = 12;
  class_count[30] = 10;

  EXPECT_EQ(10.1, mostFrequentValue(class_count, random_number_generator));
}

TEST(mostFrequentValue, equal1) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::unordered_map<double, size_t> class_count;
  class_count[1] = 10;
  class_count[2] = 15;
  class_count[3] = 15;

  EXPECT_LE(2, mostFrequentValue(class_count, random_number_generator));
  EXPECT_GE(3, mostFrequentValue(class_count, random_number_generator));
}

TEST(mostFrequentValue, equal2) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::unordered_map<double, size_t> class_count;
  class_count[10] = 30;
  class_count[11] = 30;
  class_count[15] = 29;

  EXPECT_LE(10, mostFrequentValue(class_count, random_number_generator));
  EXPECT_GE(11, mostFrequentValue(class_count, random_number_generator));
}

TEST(mostFrequentValue, equal3) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::unordered_map<double, size_t> class_count;
  class_count[3] = 10;
  class_count[5] = 500;
  class_count[6] = 500;

  EXPECT_LE(5, mostFrequentValue(class_count, random_number_generator));
  EXPECT_GE(6, mostFrequentValue(class_count, random_number_generator));
}

TEST(beautifyTime, seconds1) {
  EXPECT_EQ("0 seconds", beautifyTime(0));
}

TEST(beautifyTime, seconds2) {
  EXPECT_EQ("30 seconds", beautifyTime(30));
}

TEST(beautifyTime, minutes1) {
  EXPECT_EQ("1 minute, 0 seconds", beautifyTime(60));
}

TEST(beautifyTime, minutes2) {
  EXPECT_EQ("38 minutes, 37 seconds", beautifyTime(2317));
}

TEST(beautifyTime, hours1) {
  EXPECT_EQ("1 hour, 0 minutes, 0 seconds", beautifyTime(3600));
}

TEST(beautifyTime, hours2) {
  EXPECT_EQ("3 hours, 44 minutes, 58 seconds", beautifyTime(13498));
}

TEST(beautifyTime, days1) {
  EXPECT_EQ("1 day, 0 hours, 0 minutes, 0 seconds", beautifyTime(86400));
}

TEST(beautifyTime, days2) {
  EXPECT_EQ("3 days, 7 hours, 49 minutes, 5 seconds", beautifyTime(287345));
}

TEST(roundToNextMultiple, test0) {
  EXPECT_EQ(0, roundToNextMultiple(0, 4));
}

TEST(roundToNextMultiple, test1) {
  EXPECT_EQ(4, roundToNextMultiple(1, 4));
}

TEST(roundToNextMultiple, test2) {
  EXPECT_EQ(4, roundToNextMultiple(2, 4));
}

TEST(roundToNextMultiple, test3) {
  EXPECT_EQ(4, roundToNextMultiple(3, 4));
}

TEST(roundToNextMultiple, test4) {
  EXPECT_EQ(4, roundToNextMultiple(4, 4));
}

TEST(roundToNextMultiple, test5) {
  EXPECT_EQ(8, roundToNextMultiple(5, 4));
}

TEST(roundToNextMultiple, test6) {
  EXPECT_EQ(8, roundToNextMultiple(6, 4));
}

TEST(roundToNextMultiple, test7) {
  EXPECT_EQ(8, roundToNextMultiple(7, 4));
}

TEST(roundToNextMultiple, test8) {
  EXPECT_EQ(8, roundToNextMultiple(8, 4));
}

TEST(roundToNextMultiple, test9) {
  EXPECT_EQ(12, roundToNextMultiple(9, 4));
}

TEST(splitString, test1) {
  std::string test_string = "abc,def,ghi";
  std::vector<std::string> splitted_string;
  splitString(splitted_string, test_string, ',');
  std::vector<std::string> expect = std::vector<std::string>( { "abc", "def", "ghi" });

  EXPECT_EQ(expect, splitted_string);
}

TEST(splitString, test2) {
  std::string test_string = "abc";
  std::vector<std::string> splitted_string;
  splitString(splitted_string, test_string, ',');
  std::vector<std::string> expect = std::vector<std::string>( { "abc" });

  EXPECT_EQ(expect, splitted_string);
}

TEST(splitString, test3) {
  std::string test_string = "a-b-c";
  std::vector<std::string> splitted_string;
  splitString(splitted_string, test_string, '-');
  std::vector<std::string> expect = std::vector<std::string>( { "a", "b", "c" });

  EXPECT_EQ(expect, splitted_string);
}

TEST(splitString, test4) {
  std::string test_string = "";
  std::vector<std::string> splitted_string;
  splitString(splitted_string, test_string, ',');
  std::vector<std::string> expect = std::vector<std::string>();

  EXPECT_EQ(expect, splitted_string);
}

TEST(shuffleAndSplit, test1) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<size_t> first_part;
  std::vector<size_t> second_part;

  shuffleAndSplit(first_part, second_part, 10, 3, random_number_generator);

  EXPECT_EQ(3, first_part.size());
  EXPECT_EQ(7, second_part.size());
}

TEST(shuffleAndSplit, test2) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<size_t> first_part;
  std::vector<size_t> second_part;

  shuffleAndSplit(first_part, second_part, 100, 63, random_number_generator);

  EXPECT_EQ(63, first_part.size());
  EXPECT_EQ(37, second_part.size());
}

TEST(shuffleAndSplit, test3) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<size_t> first_part;
  std::vector<size_t> second_part;

  shuffleAndSplit(first_part, second_part, 1, 1, random_number_generator);

  EXPECT_EQ(1, first_part.size());
  EXPECT_EQ(0, second_part.size());
}

TEST(shuffleAndSplit, test4) {

  std::mt19937_64 random_number_generator;
  std::random_device random_device;
  random_number_generator.seed(random_device());

  std::vector<size_t> first_part;
  std::vector<size_t> second_part;

  shuffleAndSplit(first_part, second_part, 3, 0, random_number_generator);

  EXPECT_EQ(0, first_part.size());
  EXPECT_EQ(3, second_part.size());
}

TEST(maxstatPValueLau92, test1) {

  // From R call dput(sapply(seq(0.5, 10, by = 0.5), maxstat::pLausen92, minprop = 0.1, maxprop = 0.9))
  const std::vector<double> p_expect = std::vector<double>( { 1.0, 0.967882898076573, 0.819678995766699,
      0.463872768757117, 0.189802453892004, 0.0578438845691903, 0.0133240079314344, 0.00233924318507284,
      0.00031467775847682, 3.25492795226314e-05, 2.59527010078785e-06, 1.59801511710768e-07, 7.6090999589879e-09,
      2.80479710245055e-10, 8.01032048074225e-12, 1.77366479130538e-13, 3.04652951223938e-15, 4.06114941874027e-17,
      4.20307813816918e-19, 3.37831711514353e-21 });

  // Create sequence 0.5..10
  std::vector<double> test_b(20);
  double x = 0;
  std::generate(test_b.begin(), test_b.end(), [&]() {return x += 0.5;});

  // Compute approximation
  double minprop = 0.1;
  std::vector<double> p;
  for (auto& x : test_b) {
    p.push_back(maxstatPValueLau92(x, minprop, 1 - minprop));
  }

  // Compare with expectation
  for (size_t i = 0; i < p.size(); ++i) {
    EXPECT_NEAR(p[i], p_expect[i], p[i] * 0.05);
  }
}

