#include "TreeOptions.h"


TreeOptions::TreeOptions(uint mtry,
                         uint min_node_size,
                         std::vector<double> split_select_weights,
                         std::vector<size_t> split_select_varIDs,
                         std::vector<size_t> deterministic_varIDs,
                         std::vector<size_t> no_split_variables):
  mtry(mtry),
  min_node_size(min_node_size),
  split_select_weights(split_select_weights),
  split_select_varIDs(split_select_varIDs),
  deterministic_varIDs(deterministic_varIDs),
  no_split_variables(no_split_variables) {}

const uint TreeOptions::get_mtry() const {
  return mtry;
}

const uint TreeOptions::get_min_node_size() const {
  return min_node_size;
}

std::vector<double> TreeOptions::get_split_select_weights() const {
  return split_select_weights;
}

std::vector<size_t> TreeOptions::get_split_select_varIDs() const {
  return split_select_varIDs;
}

std::vector<size_t> TreeOptions::get_deterministic_varIDs() const {
  return deterministic_varIDs;
}

std::vector<size_t> TreeOptions::get_no_split_variables() const {
  return no_split_variables;
}