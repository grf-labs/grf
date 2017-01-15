
#ifndef GRADIENTFOREST_FORESTPREDICTOR_H
#define GRADIENTFOREST_FORESTPREDICTOR_H

#include "RelabelingStrategy.h"
#include "SplittingRule.h"
#include "PredictionStrategy.h"

#include "Tree.h"
#include "TreeTrainer.h"
#include "Forest.h"

#include <thread>
#include <future>

class ForestPredictor {
public:
  ForestPredictor(uint num_threads,
                  std::shared_ptr<PredictionStrategy> prediction_strategy);

  std::vector<std::vector<double>> predict(const Forest& forest, Data* prediction_data);
  std::vector<std::vector<double>> predict_oob(const Forest& forest, Data* original_data);

private:
  std::map<size_t, std::vector<size_t>> determine_terminal_node_IDs(
      const Forest& forest,
      Data* data,
      bool oob_prediction);

  void predictTreesInThread(uint thread_idx,
                            const Forest& forest,
                            const Data *prediction_data,
                            bool oob_prediction,
                            std::promise<std::map<size_t, std::vector<size_t>>> promise);

  void add_prediction_values(size_t nodeID,
                             const PredictionValues &prediction_values,
                             std::map<std::string, double>& average_prediction_values);

  void add_sample_weights(size_t nodeID,
                          std::shared_ptr<Tree> tree,
                          std::unordered_map<size_t, double> &weights_by_sampleID);

  void normalize_prediction_values(size_t num_trees, std::map<std::string, double>& average_prediction_values);
  void normalize_sample_weights(std::unordered_map<size_t, double>& weights_by_sampleID);

  uint num_threads;
  std::vector<uint> thread_ranges;

  std::shared_ptr<PredictionStrategy> prediction_strategy;
};


#endif //GRADIENTFOREST_FORESTPREDICTOR_H
