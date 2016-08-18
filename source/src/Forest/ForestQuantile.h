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
#include <unordered_map>

#include "globals.h"
#include "Forest.h"

class ForestQuantile: public Forest {
public:
  ForestQuantile(std::vector<double>* quantiles);
  virtual ~ForestQuantile();

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
  void addSampleWeights(size_t test_sampleID,
                        size_t tree_idx,
                        std::unordered_map<size_t, double> &weights_by_sampleID);
  void normalizeSampleWeights(std::unordered_map<size_t, double> &weights_by_sampleID);
  std::vector<double> calculateQuantileCutoffs(std::unordered_map<size_t, double> &weights_by_sampleID);

  std::vector<double>* quantiles;

  DISALLOW_COPY_AND_ASSIGN(ForestQuantile);
};


#endif //RANGER_FORESTQUANTILE_H
