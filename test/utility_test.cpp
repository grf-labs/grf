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
    EXPECT_NEAR(p[i], p_expect[i], fabs(p[i] * 0.05));
  }
}

TEST(maxstatPValueLau94, test1) {

  // From R call:
//  set.seed(123)
//  N <- 50
//  m <- which(!duplicated(sort(sample(seq(0.5,10,0.5), N, replace = TRUE)))) - 1
//  dput(m)
//  dput(sapply(seq(0.5,10,0.5), maxstat::pLausen94, N = N, minprop = 0.1, maxprop = 0.1, m = m))
  std::vector<size_t> m = std::vector<size_t>( { 0, 3, 7, 8, 12, 15, 17, 18, 21, 25, 27, 30, 31, 35, 36, 39, 44, 46 });
  const std::vector<double> p_expect = std::vector<double>( { 3.24640516569147, 2.10411448384791, 1.07190625979408,
      0.426066236447909, 0.131558039703021, 0.0314590549699589, 0.00581093458428213, 0.000826972741261553,
      9.03968770946711e-05, 7.55926672751076e-06, 4.80774063093186e-07, 2.30447718702542e-08, 8.19448285148733e-10,
      2.09519635709089e-11, 3.56382497736666e-13, 3.02490041849519e-15, -4.80261133649249e-17, -1.6110288577566e-18,
      -2.62248821317204e-20, -2.84175352170915e-22 });

  // Create sequence 0.5..10
  std::vector<double> test_b(20);
  double x = 0;
  std::generate(test_b.begin(), test_b.end(), [&]() {return x += 0.5;});

  // Compute approximation
  double minprop = 0.1;
  size_t N = 50;
  std::vector<double> p;
  for (auto& x : test_b) {
    p.push_back(maxstatPValueLau94(x, minprop, 1 - minprop, N, m));
  }

  // Compare with expectation
  for (size_t i = 0; i < p.size(); ++i) {
    EXPECT_NEAR(p[i], p_expect[i], fabs(p[i] * 0.05));
  }
}

TEST(dstdnorm, test1) {

  // From R call:
  // dput(dnorm(seq(-4, 4, by = 0.5)))

  const std::vector<double> expect = std::vector<double>( { 0.000133830225764885, 0.00087268269504576,
      0.00443184841193801, 0.0175283004935685, 0.0539909665131881, 0.129517595665892, 0.241970724519143,
      0.3520653267643, 0.398942280401433, 0.3520653267643, 0.241970724519143, 0.129517595665892, 0.0539909665131881,
      0.0175283004935685, 0.00443184841193801, 0.00087268269504576, 0.000133830225764885 });

  // Create sequence -4, 4, by=0.5
  std::vector<double> test_x(17);
  double x = -4.5;
  std::generate(test_x.begin(), test_x.end(), [&]() {return x += 0.5;});

  // Compute density
  std::vector<double> density;
  for (auto& x : test_x) {
    density.push_back(dstdnorm(x));
  }

  // Compare with expectation
  for (size_t i = 0; i < density.size(); ++i) {
    EXPECT_NEAR(density[i], expect[i], fabs(density[i] * 0.05));
  }
}

TEST(pstdnorm, test1) {

  // From R call:
  // dput(pnorm(seq(-4, 4, by = 0.5)))

  const std::vector<double> expect = std::vector<double>( { 3.16712418331199e-05, 0.000232629079035525,
      0.00134989803163009, 0.00620966532577613, 0.0227501319481792, 0.0668072012688581, 0.158655253931457,
      0.308537538725987, 0.5, 0.691462461274013, 0.841344746068543, 0.933192798731142, 0.977249868051821,
      0.993790334674224, 0.99865010196837, 0.999767370920964, 0.999968328758167 });

  // Create sequence -4, 4, by=0.5
  std::vector<double> test_x(17);
  double x = -4.5;
  std::generate(test_x.begin(), test_x.end(), [&]() {return x += 0.5;});

  // Compute distribution
  std::vector<double> dist;
  for (auto& x : test_x) {
    dist.push_back(pstdnorm(x));
  }

  // Compare with expectation
  for (size_t i = 0; i < dist.size(); ++i) {
    EXPECT_NEAR(dist[i], expect[i], fabs(dist[i] * 0.05));
  }
}

TEST(adjust_pvalues, test1) {

  // From R call:
//  set.seed(123)
//  x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
//  p <- 2*pnorm(-abs(x))
//  dput(p)
//  dput(p.adjust(p, method = "BH"))

  std::vector<double> p = std::vector<double>( { 0.575155046407955, 0.817953853056094, 0.119065433061771,
      0.943789021783783, 0.897129975905795, 0.086333312789587, 0.644858724286655, 0.205849377808347, 0.492175460488491,
      0.655841439284439, 0.220921371984431, 0.718986362304244, 0.688588400057994, 0.911867952434092, 0.578319462345744,
      0.0739515289321684, 0.618589434095213, 0.0492273640304203, 0.483080935233713, 0.636362004751351,
      0.285600042559897, 0.827448656690798, 0.304889487308354, 0.466068200259841, 0.531945286062773, 0.189079625396729,
      0.000124148072892801, 0.00161395372907077, 0.0626223946736039, 2.10159024072429e-05, 0.000611494315340942,
      0.0068319089314331, 9.81478332360736e-05, 0.000105260958830543, 0.000132598801923585, 0.000225455732762676,
      0.0003795380309271, 0.00330242962313738, 0.00705922565893261, 0.00880512860101142, 0.0211501681388388,
      0.00523699658057458, 0.0828110328240387, 2.35405350797934e-07, 2.57684187404011e-05, 0.0605329775053071,
      0.00940103983967277, 0.0112979808322889, 0.000156850331255678, 0.00353834177611007 });

  const std::vector<double> expect = std::vector<double>( { 0.741435208135569, 0.880264528394466, 0.212616844753163,
      0.943789021783783, 0.930477502483768, 0.159876505165902, 0.762606324749347, 0.343082296347245, 0.683577028456238,
      0.762606324749347, 0.356324793523276, 0.798873735893604, 0.782486818247721, 0.930477502483768, 0.741435208135569,
      0.147903057864337, 0.762606324749347, 0.11188037279641, 0.683577028456238, 0.762606324749347, 0.446250066499839,
      0.880264528394466, 0.461953768649021, 0.683577028456238, 0.718844981165909, 0.325999354132291,
      0.000947134299454175, 0.00672480720446156, 0.130463322236675, 0.000429473645673352, 0.0027795196151861,
      0.0207624284086253, 0.000947134299454175, 0.000947134299454175, 0.000947134299454175, 0.00125253184868153,
      0.0018976901546355, 0.0126369349146788, 0.0207624284086253, 0.0244586905583651, 0.0503575431877115,
      0.0174566552685819, 0.159251986200075, 1.17702675398967e-05, 0.000429473645673352, 0.130463322236675,
      0.0247395785254547, 0.0282449520807222, 0.000980314570347985, 0.0126369349146788 });

  // Adjust p-values
  std::vector<double> adjusted = adjust_pvalues(p);

  // Compare with expectation
  for (size_t i = 0; i < p.size(); ++i) {
    EXPECT_NEAR(adjusted[i], expect[i], fabs(adjusted[i] * 0.05));
  }
}

TEST(order, test1) {

  // From R call:
  //  set.seed(123)
  //  x <- runif(50)
  //  dput(x)
  //  dput(order(x, decreasing = FALSE) - 1, control = NULL)
  //  dput(order(x, decreasing = TRUE) - 1, control = NULL)

  std::vector<double> x = std::vector<double>( { 0.287577520124614, 0.788305135443807, 0.4089769218117,
      0.883017404004931, 0.940467284293845, 0.0455564993899316, 0.528105488047004, 0.892419044394046, 0.551435014465824,
      0.456614735303447, 0.956833345349878, 0.453334156190977, 0.677570635452867, 0.572633401956409, 0.102924682665616,
      0.899824970401824, 0.24608773435466, 0.0420595335308462, 0.327920719282702, 0.954503649147227, 0.889539316063747,
      0.6928034061566, 0.640506813768297, 0.994269776623696, 0.655705799115822, 0.708530468167737, 0.544066024711356,
      0.59414202044718, 0.28915973729454, 0.147113647311926, 0.963024232536554, 0.902299045119435, 0.690705278422683,
      0.795467417687178, 0.0246136845089495, 0.477795971091837, 0.758459537522867, 0.216407935833558, 0.318181007634848,
      0.231625785352662, 0.142800022382289, 0.414546335814521, 0.413724326295778, 0.368845450924709, 0.152444747742265,
      0.13880606344901, 0.233034099452198, 0.465962450252846, 0.265972640365362, 0.857827715342864 });

  const std::vector<double> expect_inc = std::vector<double>( { 34, 17, 5, 14, 45, 40, 29, 44, 37, 39, 46, 16, 48, 0,
      28, 38, 18, 43, 2, 42, 41, 11, 9, 47, 35, 6, 26, 8, 13, 27, 22, 24, 12, 32, 21, 25, 36, 1, 33, 49, 3, 20, 7, 15,
      31, 4, 19, 10, 30, 23 });

  const std::vector<double> expect_dec = std::vector<double>( { 23, 30, 10, 19, 4, 31, 15, 7, 20, 3, 49, 33, 1, 36, 25,
      21, 32, 12, 24, 22, 27, 13, 8, 26, 6, 35, 47, 9, 11, 41, 42, 2, 43, 18, 38, 28, 0, 48, 16, 46, 39, 37, 44, 29, 40,
      45, 14, 5, 17, 34 });

  // Order
  std::vector<size_t> inc = order(x, false);
  std::vector<size_t> dec = order(x, true);

  // Compare with expectation
  for (size_t i = 0; i < x.size(); ++i) {
    EXPECT_NEAR(inc[i], expect_inc[i], fabs(inc[i] * 0.05));
    EXPECT_NEAR(dec[i], expect_dec[i], fabs(dec[i] * 0.05));
  }
}

