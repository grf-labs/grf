#include "TreeOptions.h"


TreeOptions::TreeOptions(uint mtry,
                         uint min_node_size,
                         std::string split_select_weights_file,
                         std::vector<std::string> &always_split_variable_names):
  mtry(mtry),
  min_node_size(min_node_size),
  split_select_weights_file(split_select_weights_file),
  always_split_variable_names(always_split_variable_names) {}

const uint TreeOptions::get_mtry() const {
  return mtry;
}

const uint TreeOptions::get_min_node_size() const {
  return min_node_size;
}

const std::string TreeOptions::get_split_select_weights_file() const {
  return split_select_weights_file;
}

const std::vector<std::string> TreeOptions::get_always_split_variable_names() const {
  return always_split_variable_names;
}
