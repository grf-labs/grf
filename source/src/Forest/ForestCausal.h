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
  std::vector<std::vector<double>> quantile_predictions;

  DISALLOW_COPY_AND_ASSIGN(ForestCausal);
};


#endif //RANGER_FORESTCAUSAL_H
