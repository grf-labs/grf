#include <string>
#include "RegressionPredictionStrategy.h"


size_t RegressionPredictionStrategy::prediction_length() {
    return 1;
}

std::vector<double> RegressionPredictionStrategy::predict(const std::map<std::string, double>& average_prediction_values,
                                                          const std::unordered_map<size_t, double>& weights_by_sampleID,
                                                          const Observations& observations) {
  return { average_prediction_values.at(PredictionValues::AVERAGE) };
}

bool RegressionPredictionStrategy::requires_leaf_sampleIDs() {
  return false;
}

PredictionValues RegressionPredictionStrategy::precompute_prediction_values(
    const std::vector<std::vector<size_t>>& leaf_sampleIDs,
    const Observations& observations) {

  std::vector<double> averages;
  for (auto& leaf_node : leaf_sampleIDs) {
    if (leaf_node.empty()) {
      averages.push_back(0.0);
      continue;
    }

    double average = 0.0;
    for (auto& sampleID : leaf_node) {
      average += observations.get(Observations::OUTCOME).at(sampleID);
    }
    averages.push_back(average / leaf_node.size());
  }

  return PredictionValues({{PredictionValues::AVERAGE, averages}}, leaf_sampleIDs.size());
}
