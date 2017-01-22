#include <cmath>
#include <string>
#include "RegressionPredictionStrategy.h"

size_t RegressionPredictionStrategy::prediction_length() {
    return 1;
}

Prediction RegressionPredictionStrategy::predict(const std::map<std::string, double>& average_prediction_values,
                                                 const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                 const Observations& observations) {
  std::vector<double> predictions = { average_prediction_values.at(PredictionValues::AVERAGE) };
  return Prediction(predictions);
}

Prediction RegressionPredictionStrategy::predict_with_variance(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations,
    uint ci_group_size) {
  throw std::runtime_error("Variance estimates are not yet implemented.");
}

bool RegressionPredictionStrategy::requires_leaf_sampleIDs() {
  return false;
}

PredictionValues RegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations) {

  std::map<std::string, std::vector<double>> values;
  std::vector<double>& averages = values[PredictionValues::AVERAGE];

  for (auto& leaf_node : leaf_sampleIDs) {
    if (leaf_node.empty()) {
      averages.push_back(NAN);
      continue;
    }

    double average = 0.0;
    for (auto& sampleID : leaf_node) {
      average += observations.get(Observations::OUTCOME).at(sampleID);
    }
    averages.push_back(average / leaf_node.size());
  }

  return PredictionValues(values, leaf_sampleIDs.size());
}
