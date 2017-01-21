#include "ForestOptions.h"

ForestOptions::ForestOptions(uint num_trees, uint num_threads, uint random_seed):
    num_trees(num_trees), num_threads(num_threads), random_seed(random_seed) {}

uint ForestOptions::get_num_trees() {
  return num_trees;
}

uint ForestOptions::get_num_threads() {
  return num_threads;
}

uint ForestOptions::get_random_seed() {
  return random_seed;
}
