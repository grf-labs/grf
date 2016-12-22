#ifndef GRADIENTFOREST_OBSERVATIONS_H
#define GRADIENTFOREST_OBSERVATIONS_H


#include <vector>
#include <string>
#include <unordered_map>

class Observations {
public:
  Observations(std::unordered_map<std::string, std::vector<double>> observationsByType,
               size_t size);
  std::vector<double> get(std::string type);

  const size_t get_num_samples() const {
    return num_samples;
  }

  static const std::string OUTCOME = "outcome";
  static const std::string TREATMENT = "treatment";
  static const std::string INSTRUMENT = "instrument";

private:
  std::unordered_map<std::string, std::vector<double>> observationsByType;
  size_t num_samples;
};


#endif //GRADIENTFOREST_OBSERVATIONS_H
