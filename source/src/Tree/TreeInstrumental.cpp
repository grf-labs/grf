
#include <set>
#include <unordered_map>
#include "utility.h"
#include "TreeInstrumental.h"

TreeInstrumental::TreeInstrumental(size_t treatment_varID, size_t instrument_varID, std::string instrument_var_name) :
    treatment_varID(treatment_varID),
    instrument_varID(instrument_varID),
    instrument_var_name(instrument_var_name),
    udist(std::uniform_int_distribution<uint>()) {}

TreeInstrumental::TreeInstrumental(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
                                   std::vector<double> &split_values,
                                   std::vector<std::vector<size_t>> sampleIDs,
                                   size_t treatment_varID, size_t instrument_varID,
                                   std::string instrument_var_name) :
    TreeRegression(child_nodeIDs, split_varIDs, split_values),
    treatment_varID(treatment_varID),
    instrument_varID(instrument_varID),
    instrument_var_name(instrument_var_name),
    udist(std::uniform_int_distribution<uint>()) {
  this->sampleIDs = sampleIDs;
}

bool TreeInstrumental::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {

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
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  std::unordered_map<size_t, double> responses_by_sampleIDs = relabelResponses(sampleIDs[nodeID]);

  // This is an ugly hack. Should not allow splits on treatment var, etc.
  std::vector<size_t> actually_possible_split_varIDs;
  for (size_t candidate_varID : possible_split_varIDs) {
    if (candidate_varID != instrument_varID && candidate_varID != treatment_varID && candidate_varID != dependent_varID && candidate_varID != 0) {
      actually_possible_split_varIDs.push_back(candidate_varID);
    }
  }

  bool stop = responses_by_sampleIDs.empty() || findBestSplit(nodeID, actually_possible_split_varIDs, responses_by_sampleIDs);

  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }
  return false;
}

std::vector<size_t> TreeInstrumental::get_neighboring_samples(size_t sampleID) {
  size_t nodeID = prediction_terminal_nodeIDs[sampleID];
  return sampleIDs[nodeID];
}

std::unordered_map<size_t, double> TreeInstrumental::relabelResponses(std::vector<size_t>& node_sampleIDs) {
  // Calculate the relevant averages;
  size_t num_samples = node_sampleIDs.size();

  double total_treatment = 0.0;
  double total_instrument = 0.0;
  double total_response = 0.0;

  for (size_t sampleID : node_sampleIDs) {
    total_treatment += data->get(sampleID, treatment_varID);
    total_instrument += data->get(sampleID, instrument_varID);
    total_response += data->get(sampleID, dependent_varID);
  }

  double average_treatment = total_treatment / num_samples;
  double average_instrument = total_instrument / num_samples;
  double average_response = total_response / num_samples;

  // Calculate the treatment effect.
  double numerator = 0.0;
  double denominator = 0.0;
  double regularizer = 0.0;
  for (size_t sampleID : node_sampleIDs) {
    double treatment = data->get(sampleID, treatment_varID);
    double instrument = data->get(sampleID, instrument_varID);
    double response = data->get(sampleID, dependent_varID);

    numerator += (instrument - average_instrument) * (response - average_response);
    denominator += (instrument - average_instrument) * (treatment - average_treatment);
    regularizer += (instrument - average_instrument) * (instrument - average_instrument);
  }

  if (equalDoubles(denominator, 0.0)) {
    return std::unordered_map<size_t, double>(); // Signals that we should not perform a split.
  }

  // 50 is a hack
  double local_average_treatment_effect = numerator / (denominator + 50 / node_sampleIDs.size() * regularizer * sgn(denominator));

  // Create the new responses;
  std::unordered_map<size_t, double> new_responses_by_sampleID;

  for (size_t sampleID : node_sampleIDs) {
    double treatment = data->get(sampleID, treatment_varID);
    double instrument = data->get(sampleID, instrument_varID);
    double response = data->get(sampleID, dependent_varID);

    double residual = (response - average_response) - local_average_treatment_effect * (treatment - average_treatment);
    new_responses_by_sampleID[sampleID] = (instrument - average_instrument) * residual;
  }
  return new_responses_by_sampleID;
}

bool TreeInstrumental::equalDoubles(double first, double second) {
  return std::abs(first - second) < std::numeric_limits<double>::epsilon();
}

void TreeInstrumental::appendToFileInternal(std::ofstream& file) {
  saveVector2D(sampleIDs, file);
}

int TreeInstrumental::sgn(double val) {
    return (0 < val) - (val < 0);
}
