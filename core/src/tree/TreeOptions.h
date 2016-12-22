#ifndef GRADIENTFOREST_TREEOPTIONS_H
#define GRADIENTFOREST_TREEOPTIONS_H


#include <globals.h>
#include <string>
#include <vector>

class TreeOptions {
public:
  TreeOptions(uint mtry,
              uint min_node_size,
              std::string split_select_weights_file,
              std::vector<std::string>& always_split_variable_names);

  const uint get_mtry() const;
  const uint get_min_node_size() const;
  const std::string get_split_select_weights_file() const;
  const std::vector<std::string> get_always_split_variable_names() const;

private:
  uint mtry;
  uint min_node_size;
  std::string split_select_weights_file;
  std::vector<std::string>& always_split_variable_names;
};

#endif //GRADIENTFOREST_TREEOPTIONS_H
