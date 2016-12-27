#ifndef GRADIENTFOREST_OBSERVATIONS_H
#define GRADIENTFOREST_OBSERVATIONS_H


#include <vector>
#include <string>
#include <unordered_map>

class Observations {
public:
  Observations(std::unordered_map<std::string, std::vector<double>> observations_by_type,
               size_t num_samples);
  std::vector<double> get(std::string type);

  const size_t get_num_samples() const {
    return num_samples;
  }

  const std::unordered_map<std::string, std::vector<double>> get_observations_by_type() const {
    return observations_by_type;
  }

  static const std::string OUTCOME;
  static const std::string TREATMENT;
  static const std::string INSTRUMENT;

private:
  std::unordered_map<std::string, std::vector<double>> observations_by_type;
  size_t num_samples;
};


#endif //GRADIENTFOREST_OBSERVATIONS_H
