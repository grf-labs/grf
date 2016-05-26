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

TEST(adjustPvalues, test1) {

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
  std::vector<double> adjusted = adjustPvalues(p);

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

TEST(logrankScores, test1) {

  // From R call:
  //  library(survival)
  //  library(exactRankTests)
  //  y <- Surv(veteran$time, veteran$status)
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(exactRankTests::cscores(y), control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  const std::vector<double> expect = std::vector<double>( { 0.338662090461661, -2.02010051669164, -0.705694192881291,
      -0.108471626927477, -0.0608255690715493, 0.891953239598987, 0.280039238335736, 0.0692894105115154,
      -1.20022317181429, -0.87187197712046, 0.578363772221272, 0.908082271857051, -0.26880467541842, -0.315498831477494,
      0.883823158298174, 0.64336274345551, -1.58557670716783, 0.962905226680006, 0.407180963376637, 0.850487509989263,
      -1.08408138302504, -0.817639172887656, -0.328519470427333, 0.38033411652979, -0.0380982963442766,
      0.824772803240074, -0.298216440124302, 0.733731937753275, 0.393847630043303, 0.743165900017426, 0.798224130673702,
      -0.185428442042187, 0.761857488802473, 0.622086147710829, 0.446317265670834, -1.12330009489122, 0.798224130673702,
      0.484291949215137, -0.0840813830250378, 0.674400158421496, 0.407180963376637, 0.939828303603083,
      0.352746597503915, -1.85343385002497, 0.891953239598987, 0.908082271857051, 0.216259132197089, 0.600463894321395,
      -0.0380982963442766, -0.133471626927477, 0.867294232678339, -0.426436137093999, 0.970424023672487,
      0.182360827112344, -0.496633181428975, -0.426436137093999, -0.660239647426745, -2.47010051669164,
      -0.985204856795979, 0.867294232678339, -0.922704856795979, -0.575094719890513, -0.359769470427333,
      -1.49663318142897, -0.240233246846991, 0.0892894105115153, 0.10889725364877, -0.863881327384215, 0.12812802287954,
      -4.30343385002497, 0.00634614810016787, -0.751218392921895, -1.75569419288129, -0.808325771828659,
      -3.30343385002497, 0.0276227438448486, 0.985294117647059, -2.80343385002497, -1.71057670716783, 0.611333459538786,
      0.684501168522506, -1.37446559605672, -2.22010051669164, -0.61676138655718, 0.985294117647059, 0.64336274345551,
      0.555764180801456, -1.05187152346265, 0.833393492895246, 0.684501168522506, -0.89110274635123, 0.743165900017426,
      0.850487509989263, 0.248781607078105, 0.977886710239651, 0.761857488802473, 0.939828303603083, 0.714501168522506,
      0.146646541398058, 0.908082271857051, 0.146646541398058, 0.366635486392803, 0.684501168522506, 0.182360827112344,
      0.294964611470064, 0.484291949215137, 0.664196076788843, 0.714501168522506, 0.798224130673702, -0.719960761664264,
      0.622086147710829, 0.484291949215137, 0.23265257482004, 0.446317265670834, 0.324376376175947, 0.908082271857051,
      0.589474883332384, 0.532642020951149, 0.939828303603083, -0.212455469069214, -0.535094719890513,
      0.264654622951121, 0.780206112655684, 0.544269927927893, 0.294964611470064, 0.446317265670834, -0.460918895714689,
      0.780206112655684, 0.433496752850321, 0.833393492895246, 0.567127817165093, -1.28355650514763, -0.159112652568502,
      0.0276227438448486, -0.75569419288129, -1.47446559605672, 0.520877315068796 });

  // Order
  std::vector<double> scores = logrankScores(time, status);

  // Compare with expectation
  for (size_t i = 0; i < time.size(); ++i) {
    EXPECT_NEAR(scores[i], expect[i], fabs(scores[i] * 0.05));
  }
}

TEST(maxstat, trt) {

  // From R call:
  //  library(survival)
  //  library(maxstat)
  //  y <- Surv(veteran$time, veteran$status)
  //  x <- veteran$trt
  //  m <- maxstat(y, x, pmethod = "Lau92", smethod = "LogRank")
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(x, control = NULL)
  //  dput(m$statistic, control = NULL)
  //  dput(m$estimate, control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  std::vector<double> x = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 });

  const double expect_maxstat = 0.095071135385996;
  const double expect_split = 1;

  // Order
  std::vector<size_t> indices = order(x, false);

  // Scores
  std::vector<double> scores = logrankScores(time, status);

  double best_maxstat;
  double best_split_value;
  maxstat(scores, x, indices, best_maxstat, best_split_value, 0.1, 0.9);

  // Compare with expectation
  EXPECT_NEAR(best_maxstat, expect_maxstat, fabs(best_maxstat * 0.05));
  EXPECT_NEAR(best_split_value, expect_split, fabs(best_split_value * 0.05));
}

TEST(maxstat, celltype) {

  // From R call:
  //  library(survival)
  //  library(maxstat)
  //  y <- Surv(veteran$time, veteran$status)
  //  x <- veteran$trt
  //  m <- maxstat(y, x, pmethod = "Lau92", smethod = "LogRank")
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(x, control = NULL)
  //  dput(m$statistic, control = NULL)
  //  dput(m$estimate, control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  std::vector<double> x = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 });

  const double expect_maxstat = 3.44263138354529;
  const double expect_split = 1;

  // Order
  std::vector<size_t> indices = order(x, false);

  // Scores
  std::vector<double> scores = logrankScores(time, status);

  double best_maxstat;
  double best_split_value;
  maxstat(scores, x, indices, best_maxstat, best_split_value, 0.1, 0.9);

  // Compare with expectation
  EXPECT_NEAR(best_maxstat, expect_maxstat, fabs(best_maxstat * 0.05));
  EXPECT_NEAR(best_split_value, expect_split, fabs(best_split_value * 0.05));
}

TEST(maxstat, karno) {

  // From R call:
  //  library(survival)
  //  library(maxstat)
  //  y <- Surv(veteran$time, veteran$status)
  //  x <- veteran$age
  //  m <- maxstat(y, x, pmethod = "Lau92", smethod = "LogRank")
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(x, control = NULL)
  //  dput(m$statistic, control = NULL)
  //  dput(m$estimate, control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  std::vector<double> x = std::vector<double>( { 60, 70, 60, 60, 70, 20, 40, 80, 50, 70, 60, 40, 30, 80, 70, 60, 60, 40,
      80, 60, 40, 60, 60, 30, 80, 30, 50, 60, 80, 40, 20, 80, 30, 75, 70, 60, 30, 60, 80, 60, 70, 50, 50, 40, 40, 20,
      70, 40, 80, 80, 50, 80, 30, 80, 50, 80, 50, 70, 60, 40, 80, 80, 70, 90, 90, 80, 80, 70, 60, 90, 80, 80, 50, 50,
      70, 70, 20, 60, 90, 30, 20, 70, 90, 80, 50, 70, 60, 90, 50, 30, 70, 20, 30, 60, 40, 30, 20, 60, 70, 80, 85, 70,
      70, 70, 50, 30, 40, 40, 40, 99, 80, 60, 60, 60, 60, 50, 70, 10, 40, 70, 90, 80, 50, 40, 40, 60, 70, 30, 60, 30,
      60, 80, 75, 60, 70, 80, 30 });

  const double expect_maxstat = 4.61806159115936;
  const double expect_split = 40;

  // Order
  std::vector<size_t> indices = order(x, false);

  // Scores
  std::vector<double> scores = logrankScores(time, status);

  double best_maxstat;
  double best_split_value;
  maxstat(scores, x, indices, best_maxstat, best_split_value, 0.1, 0.9);

  // Compare with expectation
  EXPECT_NEAR(best_maxstat, expect_maxstat, fabs(best_maxstat * 0.05));
  EXPECT_NEAR(best_split_value, expect_split, fabs(best_split_value * 0.05));
}

TEST(maxstat, diagtime) {

  // From R call:
  //  library(survival)
  //  library(maxstat)
  //  y <- Surv(veteran$time, veteran$status)
  //  x <- veteran$age
  //  m <- maxstat(y, x, pmethod = "Lau92", smethod = "LogRank")
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(x, control = NULL)
  //  dput(m$statistic, control = NULL)
  //  dput(m$estimate, control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  std::vector<double> x = std::vector<double>( { 7, 5, 3, 9, 11, 5, 10, 29, 18, 6, 4, 58, 4, 9, 11, 3, 9, 2, 4, 4, 3, 5,
      14, 2, 3, 4, 12, 4, 12, 2, 15, 2, 5, 3, 2, 25, 4, 1, 28, 8, 1, 7, 11, 4, 23, 19, 10, 6, 2, 5, 4, 5, 3, 4, 16, 5,
      15, 2, 12, 12, 5, 12, 2, 2, 8, 11, 5, 8, 13, 12, 6, 3, 8, 1, 7, 3, 21, 3, 2, 6, 36, 13, 2, 28, 7, 11, 13, 2, 13,
      2, 22, 4, 2, 2, 36, 9, 11, 8, 3, 2, 4, 2, 2, 1, 17, 87, 8, 2, 5, 3, 3, 5, 22, 3, 3, 5, 8, 4, 4, 3, 3, 4, 10, 3, 4,
      4, 15, 4, 12, 5, 11, 10, 1, 5, 18, 4, 3 });

  const double expect_maxstat = 0.800489478294775;
  const double expect_split = 3;

  // Order
  std::vector<size_t> indices = order(x, false);

  // Scores
  std::vector<double> scores = logrankScores(time, status);

  double best_maxstat;
  double best_split_value;
  maxstat(scores, x, indices, best_maxstat, best_split_value, 0.1, 0.9);

  // Compare with expectation
  EXPECT_NEAR(best_maxstat, expect_maxstat, fabs(best_maxstat * 0.05));
  EXPECT_NEAR(best_split_value, expect_split, fabs(best_split_value * 0.05));
}

TEST(maxstat, age) {

  // From R call:
  //  library(survival)
  //  library(maxstat)
  //  y <- Surv(veteran$time, veteran$status)
  //  x <- veteran$age
  //  m <- maxstat(y, x, pmethod = "Lau92", smethod = "LogRank")
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(x, control = NULL)
  //  dput(m$statistic, control = NULL)
  //  dput(m$estimate, control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  std::vector<double> x = std::vector<double>( { 69, 64, 38, 63, 65, 49, 69, 68, 43, 70, 81, 63, 63, 52, 48, 61, 42, 35,
      63, 56, 55, 67, 63, 65, 46, 53, 69, 68, 43, 55, 42, 64, 65, 65, 55, 66, 60, 67, 53, 62, 67, 72, 48, 68, 67, 61,
      60, 62, 38, 50, 63, 64, 43, 34, 66, 62, 52, 47, 63, 68, 45, 41, 66, 62, 60, 66, 38, 53, 37, 54, 60, 48, 52, 70,
      50, 62, 65, 58, 62, 64, 63, 58, 64, 52, 35, 63, 70, 51, 40, 69, 36, 71, 62, 60, 44, 54, 66, 49, 72, 68, 62, 71,
      70, 61, 71, 59, 67, 60, 69, 57, 39, 62, 50, 43, 70, 66, 61, 81, 58, 63, 60, 62, 42, 69, 63, 45, 68, 39, 66, 63,
      49, 64, 65, 64, 67, 65, 37 });

  const double expect_maxstat = 1.7992993341166;
  const double expect_split = 58;

  // Order
  std::vector<size_t> indices = order(x, false);

  // Scores
  std::vector<double> scores = logrankScores(time, status);

  double best_maxstat;
  double best_split_value;
  maxstat(scores, x, indices, best_maxstat, best_split_value, 0.1, 0.9);

  // Compare with expectation
  EXPECT_NEAR(best_maxstat, expect_maxstat, fabs(best_maxstat * 0.05));
  EXPECT_NEAR(best_split_value, expect_split, fabs(best_split_value * 0.05));
}

TEST(maxstat, prior) {

  // From R call:
  //  library(survival)
  //  library(maxstat)
  //  y <- Surv(veteran$time, veteran$status)
  //  x <- veteran$age
  //  m <- maxstat(y, x, pmethod = "Lau92", smethod = "LogRank")
  //  dput(y[, 1], control = NULL)
  //  dput(y[, 2], control = NULL)
  //  dput(x, control = NULL)
  //  dput(m$statistic, control = NULL)
  //  dput(m$estimate, control = NULL)

  std::vector<double> time = std::vector<double>( { 72, 411, 228, 126, 118, 10, 82, 110, 314, 100, 42, 8, 144, 25, 11,
      30, 384, 4, 54, 13, 123, 97, 153, 59, 117, 16, 151, 22, 56, 21, 18, 139, 20, 31, 52, 287, 18, 51, 122, 27, 54, 7,
      63, 392, 10, 8, 92, 35, 117, 132, 12, 162, 3, 95, 177, 162, 216, 553, 278, 12, 260, 200, 156, 182, 143, 105, 103,
      250, 100, 999, 112, 87, 231, 242, 991, 111, 1, 587, 389, 33, 25, 357, 467, 201, 1, 30, 44, 283, 15, 25, 103, 21,
      13, 87, 2, 20, 7, 24, 99, 8, 99, 61, 25, 95, 80, 51, 29, 24, 18, 83, 31, 51, 90, 52, 73, 8, 36, 48, 7, 140, 186,
      84, 19, 45, 80, 52, 164, 19, 53, 15, 43, 340, 133, 111, 231, 378, 49 });

  std::vector<double> status = std::vector<double>( { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1 });

  std::vector<double> x = std::vector<double>( { 0, 10, 0, 10, 10, 0, 10, 0, 0, 0, 0, 10, 0, 10, 10, 0, 0, 0, 10, 0, 0,
      0, 10, 0, 0, 10, 0, 0, 10, 10, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 0, 0, 0, 0, 10, 0, 0, 0, 10, 0,
      0, 0, 0, 10, 0, 10, 0, 0, 0, 0, 0, 10, 10, 10, 0, 0, 10, 0, 10, 0, 10, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10, 0, 10, 0,
      10, 0, 0, 0, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 10, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0,
      10, 10, 0, 0, 10, 10, 0, 0, 10, 0, 0 });

  const double expect_maxstat = 0.7158562778385;
  const double expect_split = 0;

  // Order
  std::vector<size_t> indices = order(x, false);

  // Scores
  std::vector<double> scores = logrankScores(time, status);

  double best_maxstat;
  double best_split_value;
  maxstat(scores, x, indices, best_maxstat, best_split_value, 0.1, 0.9);

  // Compare with expectation
  EXPECT_NEAR(best_maxstat, expect_maxstat, fabs(best_maxstat * 0.05));
  EXPECT_NEAR(best_split_value, expect_split, fabs(best_split_value * 0.05));
}

TEST(numSamplesLeftOfCutpoint, test1) {

  // From R call:
  //  library(survival)
  //  x <- veteran$age
  //  ties <- duplicated(sort(x))
  //  m <- (which(!ties) - 1)[-1]
  //  if (ties[length(x)]) {
  //    m <- c(m, length(x))
  //  }
  //  dput(x, control = NULL)
  //  dput(m, control = NULL)

  std::vector<double> x = std::vector<double>( { 69, 64, 38, 63, 65, 49, 69, 68, 43, 70, 81, 63, 63, 52, 48, 61, 42, 35,
      63, 56, 55, 67, 63, 65, 46, 53, 69, 68, 43, 55, 42, 64, 65, 65, 55, 66, 60, 67, 53, 62, 67, 72, 48, 68, 67, 61,
      60, 62, 38, 50, 63, 64, 43, 34, 66, 62, 52, 47, 63, 68, 45, 41, 66, 62, 60, 66, 38, 53, 37, 54, 60, 48, 52, 70,
      50, 62, 65, 58, 62, 64, 63, 58, 64, 52, 35, 63, 70, 51, 40, 69, 36, 71, 62, 60, 44, 54, 66, 49, 72, 68, 62, 71,
      70, 61, 71, 59, 67, 60, 69, 57, 39, 62, 50, 43, 70, 66, 61, 81, 58, 63, 60, 62, 42, 69, 63, 45, 68, 39, 66, 63,
      49, 64, 65, 64, 67, 65, 37 });

  std::vector<size_t> expect = std::vector<size_t>( { 1, 3, 4, 6, 9, 11, 12, 13, 16, 20, 21, 23, 24, 25, 28, 31, 34, 35,
      39, 42, 44, 47, 48, 49, 52, 53, 60, 64, 74, 86, 93, 100, 107, 113, 119, 125, 130, 133, 135, 137 });

  // Order
  std::vector<size_t> indices = order(x, false);

  std::vector<size_t> m = numSamplesLeftOfCutpoint(x, indices);

  // Compare with expectation
  EXPECT_EQ(m.size(), expect.size());
  for (size_t i = 0; i < expect.size(); ++i) {
    EXPECT_EQ(m[i], expect[i]);
  }
}

