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
#include "TreeFactory.h"

class ForestQuantile : public Forest {
public:
  ForestQuantile(std::vector<double> *quantiles,
                 std::unordered_map<std::string, size_t> observables,
                 RelabelingStrategy *relabeling_strategy,
                 SplittingRule *splitting_rule);

  virtual ~ForestQuantile();

private:
  void initInternal(std::string status_variable_name);
  void computePredictionErrorInternal();
  void writeOutputInternal();
  void writeConfusionFile();
  void writePredictionFile();
  void saveToFileInternal(std::ofstream& outfile);
  void loadFromFileInternal(std::ifstream& infile);

  void predictInternal();
  void addSampleWeights(size_t test_sample_idx,
                        TreeFactory* tree,
                        std::unordered_map<size_t, double> &weights_by_sampleID);
  void normalizeSampleWeights(std::unordered_map<size_t, double> &weights_by_sampleID);

  std::vector<double>* quantiles;

  DISALLOW_COPY_AND_ASSIGN(ForestQuantile);
};


#endif //RANGER_FORESTQUANTILE_H
