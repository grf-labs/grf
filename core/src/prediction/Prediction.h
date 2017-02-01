#ifndef GRADIENTFOREST_PREDICTIONRESULTS_H
#define GRADIENTFOREST_PREDICTIONRESULTS_H

#include <cstddef>
#include <vector>

class Prediction {
public:
  Prediction(const std::vector<double>& predictions);
  Prediction(const std::vector<double>& predictions,
             const std::vector<double>& variance_estimates);

  const std::vector<double>& get_predictions() {
    return predictions;
  }

  const std::vector<double>& get_variance_estimates() {
    return variance_estimates;
  }

  bool contains_variance_estimates() {
    return !variance_estimates.empty();
  }

  size_t size() {
    return predictions.size();
  }

private:
  std::vector<double> predictions;
  std::vector<double> variance_estimates;
};


#endif //GRADIENTFOREST_PREDICTIONRESULTS_H
