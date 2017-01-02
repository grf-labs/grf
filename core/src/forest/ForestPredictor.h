
#ifndef GRADIENTFOREST_FORESTPREDICTOR_H
#define GRADIENTFOREST_FORESTPREDICTOR_H

#include "RelabelingStrategy.h"
#include "SplittingRule.h"
#include "PredictionStrategy.h"

#include "Tree.h"
#include "TreeModel.h"
#include "Forest.h"

#include <thread>
#include <future>

class ForestPredictor {
public:
  ForestPredictor(PredictionStrategy *prediction_strategy);

  std::vector<std::vector<double>> predict(Forest* forest, Data* prediction_data);

  void init(std::string load_forest_filename,
            uint num_threads,
            std::ostream *verbose_out);

  void writeConfusionFile(Data* prediction_data, std::vector<std::vector<double>> predictions);
  void writePredictionFile(Data* prediction_data, std::vector<std::vector<double>> predictions);

private:
  void computePredictionError(Forest* forest, Data* prediction_data);
  void computePredictionErrorInternal(Forest* forest,
                                      Data* prediction_data,
                                      std::unordered_map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree);

  void predictTreesInThread(uint thread_idx,
                            Forest *forest,
                            const Data *prediction_data,
                            bool oob_prediction,
                            std::promise<std::unordered_map<size_t, std::vector<size_t>>> promise);

  void addSampleWeights(size_t test_sample_idx,
                        Tree* tree,
                        std::unordered_map<size_t, double> &weights_by_sampleID,
                        std::vector<size_t> terminal_node_IDs);
  void normalizeSampleWeights(std::unordered_map<size_t, double> &weights_by_sampleID);
  std::vector<std::vector<double>> predictInternal(Forest* forest,
                                                   Data* prediction_data,
                                                   std::unordered_map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree);

  std::ostream* verbose_out;

  bool prediction_mode;

  // Multithreading
  uint num_threads;
  std::vector<uint> thread_ranges;

  PredictionStrategy* prediction_strategy;
};


#endif //GRADIENTFOREST_FORESTPREDICTOR_H
