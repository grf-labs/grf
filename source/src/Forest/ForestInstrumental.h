#ifndef RANGER_ForestInstrumental_H
#define RANGER_ForestInstrumental_H

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
#include "InstrumentalTreeFactory.h"

class ForestInstrumental: public Forest {
public:
  ForestInstrumental(std::string instrument_variable_name);
  virtual ~ForestInstrumental();

private:
  void initInternal(std::string status_variable_name);
  void growInternal();
  void computePredictionErrorInternal();
  void writeOutputInternal();
  void writeConfusionFile();
  void writePredictionFile();
  void saveToFileInternal(std::ofstream& outfile);
  void loadFromFileInternal(std::ifstream& infile);

  void addSampleWeights(size_t test_sample_idx,
                        InstrumentalTreeFactory* tree,
                        std::unordered_map<size_t, double> &weights_by_sampleID);
  void normalizeSampleWeights(std::unordered_map<size_t, double> &weights_by_sampleID);
  void predictInternal();

  size_t treatment_varID;
  size_t instrument_varID;
  std::string instrument_variable_name;
  std::unordered_map<size_t, std::vector<double>>* original_responses;

  DISALLOW_COPY_AND_ASSIGN(ForestInstrumental);
};


#endif //RANGER_ForestInstrumental_H
