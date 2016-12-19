//
// Created by Julie Noelle Tibshirani on 12/17/16.
//

#ifndef RANGER_TREEMODEL_H
#define RANGER_TREEMODEL_H

#include "prediction/PredictionStrategy.h"
#include "relabeling/RelabelingStrategy.h"
#include "splitting/SplittingRule.h"
#include "../utility/Data.h"
#include "Tree.h"
#include "sampling/BootstrapSampler.h"

class TreeModel {
public:
  TreeModel(RelabelingStrategy *relabeling_strategy,
            SplittingRule *splitting_rule,
            PredictionStrategy *prediction_strategy,
            size_t dependent_varID,
            uint mtry,
            uint min_node_size,
            std::vector <size_t> *deterministic_varIDs,
            std::vector <size_t> *split_select_varIDs,
            std::vector <size_t> *no_split_variables);

  Tree* train(Data *data,
              std::unordered_map<std::string, std::vector<double>>* observations,
              BootstrapSampler *bootstrap_sampler,
              std::vector<double> *split_select_weights);

private:
  void createEmptyNode(std::vector <std::vector<size_t>> &child_nodeIDs,
                       std::vector <std::vector<size_t>> &sampleIDs,
                       std::vector <size_t> &split_varIDs,
                       std::vector<double> &split_values);

  void createPossibleSplitVarSubset(std::vector <size_t> &result,
                                    Data* data,
                                    std::vector<double> *split_select_weights);

  bool splitNode(size_t nodeID,
                 Data* data,
                 std::unordered_map<std::string, std::vector<double>>* observations,
                 std::vector <std::vector<size_t>> &child_nodeIDs,
                 std::vector <std::vector<size_t>> &sampleIDs,
                 std::vector <size_t> &split_varIDs,
                 std::vector<double> &split_values,
                 std::vector<double> *split_select_weights);

  bool splitNodeInternal(size_t nodeID,
                         Data* data,
                         std::unordered_map<std::string, std::vector<double>>* observations,
                         std::vector <size_t> &possible_split_varIDs,
                         std::vector <std::vector<size_t>> &sampleIDs,
                         std::vector <size_t> &split_varIDs,
                         std::vector<double> &split_values);

  RelabelingStrategy *relabeling_strategy;
  SplittingRule *splitting_rule;
  PredictionStrategy *prediction_strategy;

  size_t dependent_varID;

  uint mtry;
  uint min_node_size;

  std::vector <size_t> *deterministic_varIDs;
  std::vector <size_t> *split_select_varIDs;
  std::vector <size_t> *no_split_variables;
};

#endif //RANGER_TREEMODEL_H
