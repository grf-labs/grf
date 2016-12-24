#ifndef GRADIENTFOREST_FORESTMODEL_H
#define GRADIENTFOREST_FORESTMODEL_H

#include "RelabelingStrategy.h"
#include "SplittingRule.h"
#include "PredictionStrategy.h"

#include "Tree.h"
#include "TreeModel.h"
#include "Forest.h"

#include <thread>
#include <future>

class ForestModel {
public:
  ForestModel(std::unordered_map<std::string, size_t> observables,
              RelabelingStrategy *relabeling_strategy,
              SplittingRule *splitting_rule,
              PredictionStrategy *prediction_strategy);

  Forest* train(Data* data);
  std::vector<std::vector<double>> predict(Forest* forest, Data* prediction_data);

  void init(uint mtry,
            uint num_trees, std::ostream *verbose_out, uint seed, uint num_threads,
            std::string load_forest_filename, uint min_node_size,
            std::string split_select_weights_file, std::vector<std::string> &always_split_variable_names,
            bool sample_with_replacement,
            bool memory_saving_splitting,
            std::string case_weights_file, double sample_fraction);

  void writeConfusionFile(Data* prediction_data, std::vector<std::vector<double>> predictions);
  void writePredictionFile(Data* prediction_data, std::vector<std::vector<double>> predictions);

private:
  void computePredictionError(Forest* forest, Data* prediction_data);
  void computePredictionErrorInternal(Forest* forest,
                                      Data* prediction_data,
                                      std::unordered_map<size_t, std::vector<size_t>> terminal_node_IDs_by_tree);

  void growTreesInThread(uint thread_idx,
                         Data *data,
                         Observations* observations,
                         std::promise<std::vector<Tree*>> promise);

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

  void setSplitWeightVector(std::vector<double>& split_select_weights,
                            size_t num_independent_variables);
  void setAlwaysSplitVariables(Data* data,
                               std::vector<std::string>& always_split_variable_names,
                               size_t num_independent_variables);

  std::ostream* verbose_out;

  size_t num_trees;
  uint mtry;
  uint min_node_size;
  uint seed;
  bool prediction_mode;
  bool sample_with_replacement;
  bool memory_saving_splitting;
  bool keep_inbag;
  double sample_fraction;

  std::vector<size_t> no_split_variables;

  // Multithreading
  uint num_threads;
  std::vector<uint> thread_ranges;

  TreeModel* tree_model;

  // Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) for each variable
  // Deterministic variables are always selected
  std::vector<std::string> always_split_variable_names;
  std::vector<size_t> deterministic_varIDs;
  std::vector<size_t> split_select_varIDs;
  std::vector<double> split_select_weights;

  // Bootstrap weights
  std::vector<double> case_weights;

  // Random number generator
  std::mt19937_64 random_number_generator;

  // Computation progress (finished trees)
  size_t progress;

  std::unordered_map<std::string, size_t> observables;
  RelabelingStrategy* relabeling_strategy;
  SplittingRule* splitting_rule;
  PredictionStrategy* prediction_strategy;

  std::string split_select_weights_file;
  std::string case_weights_file;

  DISALLOW_COPY_AND_ASSIGN(ForestModel);
};



#endif //GRADIENTFOREST_FORESTMODEL_H
