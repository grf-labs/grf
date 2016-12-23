#include "TestUtilities.h"

Observations* TestUtilities::create_observations(std::vector<double> outcome) {
  std::unordered_map<std::string, std::vector<double>> observationsByType = {
      {Observations::OUTCOME, outcome}};
  return new Observations(observationsByType, outcome.size());
}

Observations* TestUtilities::create_observations(std::vector<double> outcome,
                                                 std::vector<double> treatment,
                                                 std::vector<double> instrument) {
  std::unordered_map<std::string, std::vector<double>> observationsByType = {
      {Observations::OUTCOME, outcome},
      {Observations::TREATMENT, treatment},
      {Observations::INSTRUMENT, instrument}};
  return new Observations(observationsByType, outcome.size());
}