#ifndef GRADIENTFOREST_FORESTOPTIONS_H
#define GRADIENTFOREST_FORESTOPTIONS_H


#include "globals.h"

class ForestOptions {
public:
  ForestOptions(uint num_trees, uint num_threads, uint random_seed);

  const uint get_num_trees() const;
  const uint get_num_threads() const;
  const uint get_random_seed() const;

private:
  uint num_trees;
  uint num_threads;
  uint random_seed;
};


#endif //GRADIENTFOREST_FORESTOPTIONS_H
