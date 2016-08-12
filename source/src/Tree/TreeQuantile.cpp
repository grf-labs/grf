
#include <set>
#include "TreeQuantile.h"


TreeQuantile::TreeQuantile(std::vector<double>* quantiles) : quantiles(quantiles) {}

bool TreeQuantile::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  // Check if node is pure and set split_value to estimate and stop if pure
  bool pure = true;
  double pure_value = 0;
  for (size_t i = 0; i < sampleIDs[nodeID].size(); ++i) {
    double value = data->get(sampleIDs[nodeID][i], dependent_varID);
    if (i != 0 && value != pure_value) {
      pure = false;
      break;
    }
    pure_value = value;
  }
  if (pure) {
    split_values[nodeID] = pure_value;
    return true;
  }

  std::vector<double> responses(sampleIDs[nodeID].size());
  data->getAllValues(responses, sampleIDs[nodeID], dependent_varID);
  std::vector<uint>* relabeled_responses = relabelResponses(&responses);

  TreeClassification* classificationTree = createClassificationTree(sampleIDs[nodeID],
                                                                    relabeled_responses);
  bool stop = classificationTree->findBestSplit(nodeID, possible_split_varIDs);

  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}

std::vector<uint>* TreeQuantile::relabelResponses(std::vector<double>* responses) {
  std::vector<double> sorted_responses(*responses);
  std::sort(sorted_responses.begin(), sorted_responses.end());

  size_t num_samples = responses->size();
  std::vector<double> quantile_cutoffs;

  // Calculate the response value cutoffs for each quantile.
  for (auto& quantile : *quantiles) {
    size_t response_index = (size_t) floor(num_samples * quantile);
    quantile_cutoffs.push_back(responses->at(response_index));
  }

  // Remove duplicate cutoffs.
  quantile_cutoffs.erase(std::unique(quantile_cutoffs.begin(), quantile_cutoffs.end()));

  // Assign a class to each response based on what quantile it belongs to.
  std::vector<uint>* relabeled_responses = new std::vector<uint>();
  for (auto& response : *responses) {
    auto quantile = std::lower_bound(quantile_cutoffs.begin(),
                                     quantile_cutoffs.end(),
                                     response);
    long quantile_index = quantile - quantile_cutoffs.begin();
    relabeled_responses->push_back((uint) quantile_index);
  }
  return relabeled_responses;
}

TreeClassification* TreeQuantile::createClassificationTree(std::vector<size_t>& sampleIDs,
                                                           std::vector<uint>* relabeledResponses) {
  std::set<uint>* unique_classIDs = new std::set<uint>(relabeledResponses->begin(),
                                                       relabeledResponses->end());
  size_t num_classes = unique_classIDs->size();

  // Create a dummy vector of class values. TreeClassification only uses the size of
  // this vector, so it is fine that the values are fake.
  // TODO(jtibs): Refactor TreeClassification to avoid this hack.
  std::vector<double>* class_values = new std::vector<double>(num_classes, 1.0);

  std::vector<uint>* response_classIDs = new std::vector<uint>(num_samples);
  for (size_t i = 0; i < sampleIDs.size(); ++i) {
    size_t sampleID = sampleIDs[i];
    (*response_classIDs)[sampleID] = (*relabeledResponses)[i];
  }

  TreeClassification* tree = new TreeClassification(class_values, response_classIDs);

  uint tree_seed = 42;
  tree->init(data, mtry, dependent_varID, num_samples, tree_seed, deterministic_varIDs, split_select_varIDs,
             split_select_weights, importance_mode, min_node_size, no_split_variables, sample_with_replacement,
             is_ordered_variable, memory_saving_splitting, splitrule, case_weights, keep_inbag, sample_fraction,
             alpha, minprop, holdout);
  return tree;
}