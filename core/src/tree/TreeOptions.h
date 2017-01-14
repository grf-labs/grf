#ifndef GRADIENTFOREST_TREEOPTIONS_H
#define GRADIENTFOREST_TREEOPTIONS_H


#include <string>
#include <vector>

#include "globals.h"

class TreeOptions {
public:
  TreeOptions(uint mtry,
              uint min_node_size,
              std::vector<double> split_select_weights,
              std::vector<size_t> split_select_varIDs,
              std::vector<size_t> deterministic_varIDs,
              std::vector<size_t> no_split_variables,
              bool honesty);

  const uint get_mtry() const;
  const uint get_min_node_size() const;
  std::vector<double> get_split_select_weights() const;
  std::vector<size_t> get_split_select_varIDs() const;

  std::vector<size_t> get_deterministic_varIDs() const;
  std::vector<size_t> get_no_split_variables() const;

  bool get_honesty() const;

private:
  uint mtry;
  uint min_node_size;
  std::vector<double> split_select_weights;
  std::vector<size_t> split_select_varIDs;
  std::vector<size_t> deterministic_varIDs;
  std::vector<size_t> no_split_variables;

  bool honesty;
};

#endif //GRADIENTFOREST_TREEOPTIONS_H
