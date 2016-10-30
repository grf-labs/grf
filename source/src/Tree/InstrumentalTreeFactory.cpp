
#include <set>
#include <unordered_map>
#include "utility.h"
#include "InstrumentalTreeFactory.h"
#include "InstrumentalRelabelingStrategy.h"

InstrumentalTreeFactory::InstrumentalTreeFactory(size_t treatment_varID, size_t instrument_varID, std::string instrument_var_name) :
    treatment_varID(treatment_varID),
    instrument_varID(instrument_varID),
    instrument_var_name(instrument_var_name) {}

InstrumentalTreeFactory::InstrumentalTreeFactory(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                                   std::vector<double> &split_values,
                                   std::vector<std::vector<size_t>> sampleIDs,
                                   size_t treatment_varID, size_t instrument_varID,
                                   std::string instrument_var_name) :
    TreeFactory(child_nodeIDs, split_varIDs, split_values),
    treatment_varID(treatment_varID),
    instrument_varID(instrument_varID),
    instrument_var_name(instrument_var_name) {
  this->sampleIDs = sampleIDs;
}

bool InstrumentalTreeFactory::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

  // Check node size, stop if maximum reached
  if (sampleIDs[nodeID].size() <= min_node_size) {
    split_values[nodeID] = -1.0;
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
    split_values[nodeID] = -1.0;
    return true;
  }

  InstrumentalRelabelingStrategy* relabelingStrategy = new InstrumentalRelabelingStrategy(data,
        dependent_varID, treatment_varID, instrument_varID);
  std::unordered_map<size_t, double> responses_by_sampleIDs = relabelingStrategy->relabelResponses(sampleIDs[nodeID]);

  // This is an ugly hack. Should not allow splits on treatment var, etc.
  std::vector<size_t> actually_possible_split_varIDs;
  for (size_t candidate_varID : possible_split_varIDs) {
    if (candidate_varID != instrument_varID &&
        candidate_varID != treatment_varID &&
        candidate_varID != dependent_varID &&
        candidate_varID != 0) {
      actually_possible_split_varIDs.push_back(candidate_varID);
    }
  }

  RegressionSplittingRule* splittingRule = new RegressionSplittingRule(responses_by_sampleIDs, data,
    sampleIDs, split_varIDs, split_values);
  bool stop = responses_by_sampleIDs.empty() || splittingRule->findBestSplit(nodeID, actually_possible_split_varIDs);

  if (stop) {
    split_values[nodeID] = -1.0;
    return true;
  }
  return false;
}

std::vector<size_t> InstrumentalTreeFactory::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}

