#include "Prediction.h"

Prediction::Prediction(const std::vector<double>& predictions):
  predictions(predictions), variance_estimates(0) {}

Prediction::Prediction(const std::vector<double>& predictions,
                       const std::vector<double>& variance_estimates):
  predictions(predictions), variance_estimates(variance_estimates) {}
