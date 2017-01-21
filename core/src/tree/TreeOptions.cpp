#include "TreeOptions.h"


TreeOptions::TreeOptions(uint mtry,
                         uint min_node_size,
                         const std::vector<double>& split_select_weights,
                         const std::vector<size_t>& split_select_varIDs,
                         const std::vector<size_t>& deterministic_varIDs,
                         const std::vector<size_t>& no_split_variables,
                         bool honesty):
  mtry(mtry),
  min_node_size(min_node_size),
  split_select_weights(split_select_weights),
  split_select_varIDs(split_select_varIDs),
  deterministic_varIDs(deterministic_varIDs),
  no_split_variables(no_split_variables),
  honesty(honesty) {}

uint TreeOptions::get_mtry() {
  return mtry;
}

uint TreeOptions::get_min_node_size() {
  return min_node_size;
}

const std::vector<double>& TreeOptions::get_split_select_weights() {
  return split_select_weights;
}

const std::vector<size_t>& TreeOptions::get_split_select_varIDs() {
  return split_select_varIDs;
}

const std::vector<size_t>& TreeOptions::get_deterministic_varIDs() {
  return deterministic_varIDs;
}

const std::vector<size_t>& TreeOptions::get_no_split_variables() {
  return no_split_variables;
}

bool TreeOptions::get_honesty() {
  return honesty;
}