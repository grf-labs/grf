#ifndef GRADIENTFOREST_OBSERVATIONS_H
#define GRADIENTFOREST_OBSERVATIONS_H


#include <vector>
#include <string>
#include <map>

class Observations {
public:
  Observations(const std::map<std::string, std::vector<double>>& observations_by_type,
               size_t num_samples);

  const std::vector<double>& get(std::string type) const;

  size_t get_num_samples() const;
  const std::map<std::string, std::vector<double>>& get_observations_by_type() const;

  static const std::string OUTCOME;
  static const std::string TREATMENT;
  static const std::string INSTRUMENT;

private:
  std::map<std::string, std::vector<double>> observations_by_type;
  size_t num_samples;
};


#endif //GRADIENTFOREST_OBSERVATIONS_H
