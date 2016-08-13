#ifndef RANGER_FORESTQUANTILE_H
#define RANGER_FORESTQUANTILE_H

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "globals.h"
#include "Forest.h"

#include <iostream>
#include <vector>

#include "globals.h"
#include "Forest.h"

class ForestQuantile: public Forest {
public:
  ForestQuantile(std::vector<double>* quantiles);
  virtual ~ForestQuantile();

private:
  void initInternal(std::string status_variable_name);
  void growInternal();
  void predictInternal();
  void computePredictionErrorInternal();
  void writeOutputInternal();
  void writeConfusionFile();
  void writePredictionFile();
  void saveToFileInternal(std::ofstream& outfile);
  void loadFromFileInternal(std::ifstream& infile);

  std::vector<double>* quantiles;

  DISALLOW_COPY_AND_ASSIGN(ForestQuantile);
};


#endif //RANGER_FORESTQUANTILE_H
