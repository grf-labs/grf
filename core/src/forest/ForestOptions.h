#ifndef GRADIENTFOREST_FORESTOPTIONS_H
#define GRADIENTFOREST_FORESTOPTIONS_H


#include "globals.h"

class ForestOptions {
public:
  ForestOptions(uint num_trees, uint num_threads, uint random_seed);

  uint get_num_trees();
  uint get_num_threads();
  uint get_random_seed();

private:
  uint num_trees;
  uint num_threads;
  uint random_seed;
};


#endif //GRADIENTFOREST_FORESTOPTIONS_H
