#include "ForestOptions.h"

ForestOptions::ForestOptions(uint num_trees, uint num_threads, uint random_seed):
    num_trees(num_trees), num_threads(num_threads), random_seed(random_seed) {}

const uint ForestOptions::get_num_trees() const {
  return num_trees;
}

const uint ForestOptions::get_num_threads() const {
  return num_threads;
}

const uint ForestOptions::get_random_seed() const {
  return random_seed;
}
