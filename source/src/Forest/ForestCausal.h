#ifndef RANGER_FORESTCAUSAL_H
#define RANGER_FORESTCAUSAL_H

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "globals.h"
#include "Forest.h"

#include <iostream>
#include <vector>
#include <unordered_map>

#include "globals.h"
#include "Forest.h"

class ForestCausal: public Forest {
public:
  ForestCausal();
  virtual ~ForestCausal();

  void loadForest(size_t dependent_varID, size_t num_trees,
                  std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
                  std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values);

private:
  void initInternal(std::string status_variable_name);
  void growInternal();
  void computePredictionErrorInternal();
  void writeOutputInternal();
  void writeConfusionFile();
  void writePredictionFile();
  void saveToFileInternal(std::ofstream& outfile);
  void loadFromFileInternal(std::ifstream& infile);

  void predictInternal();

  size_t treatment_varID;

  DISALLOW_COPY_AND_ASSIGN(ForestCausal);
};


#endif //RANGER_FORESTCAUSAL_H
