#include "ForestTestUtilities.h"
#include "ForestTrainer.h"

void ForestTestUtilities::init_trainer(ForestTrainer& trainer) {
  init_trainer(trainer, false);
}

void ForestTestUtilities::init_honest_trainer(ForestTrainer &trainer) {
  init_trainer(trainer, true);
}

void ForestTestUtilities::init_trainer(ForestTrainer& trainer, bool honesty) {
  uint mtry = 3;
  uint num_trees = 20;
  std::ostream* verbose_out = &std::cout;
  uint seed = 42;
  uint num_threads = 4;
  std::string load_forest_filename = "";
  uint min_node_size = 1;
  std::vector<size_t> no_split_variables;
  std::string split_select_weights_file = "";
  std::vector<std::string> always_split_variable_names;
  bool sample_with_replacement = true;
  bool memory_saving_splitting = false;
  std::string case_weights_file = "";
  double sample_fraction = 0.7;
  uint ci_bag_size = 1;

  trainer.init(mtry, num_trees, verbose_out, seed, num_threads, load_forest_filename,
               min_node_size, no_split_variables, split_select_weights_file, always_split_variable_names,
               sample_with_replacement, memory_saving_splitting, case_weights_file, sample_fraction,
               honesty, ci_bag_size);
}
