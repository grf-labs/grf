#ifndef GRADIENTFOREST_PREDICTIONSTRATEGY_H
#define GRADIENTFOREST_PREDICTIONSTRATEGY_H

#include <unordered_map>
#include <vector>

class PredictionStrategy {
public:
  virtual std::vector<double> predict(std::unordered_map<size_t, double>& weights_by_sampleID,
                                      std::unordered_map<std::string, std::vector<double>> original_observations) = 0;
};


#endif //GRADIENTFOREST_PREDICTIONSTRATEGY_H
