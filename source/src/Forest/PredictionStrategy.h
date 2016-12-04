#ifndef RANGER_PREDICTIONSTRATEGY_H
#define RANGER_PREDICTIONSTRATEGY_H


class PredictionStrategy {
public:
  virtual std::vector<double> predict(std::unordered_map<size_t, double>& weights_by_sampleID) = 0;
};


#endif //RANGER_PREDICTIONSTRATEGY_H
