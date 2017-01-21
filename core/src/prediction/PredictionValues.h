#ifndef GRADIENTFOREST_PREDICTIONVALUES_H
#define GRADIENTFOREST_PREDICTIONVALUES_H


#include <vector>
#include <string>
#include <map>

class PredictionValues {
public:
  PredictionValues();

  PredictionValues(const std::map<std::string, std::vector<double>>& values_by_type,
                   size_t num_nodes);

  const std::vector<double>& get(std::string type) const;

  const size_t get_num_nodes() const {
    return num_nodes;
  }

  const std::map<std::string, std::vector<double>>& get_values_by_type() const {
    return values_by_type;
  }

  const bool empty() const {
    return num_nodes == 0;
  }

  static const std::string AVERAGE;

private:
  std::map<std::string, std::vector<double>> values_by_type;
  size_t num_nodes;
};


#endif //GRADIENTFOREST_PREDICTIONVALUES_H
