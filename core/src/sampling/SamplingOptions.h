#ifndef GRADIENTFOREST_SAMPLINGOPTIONS_H
#define GRADIENTFOREST_SAMPLINGOPTIONS_H


#include <string>
#include <vector>

class SamplingOptions {
public:
  SamplingOptions(bool sample_with_replacement,
                  const std::vector<double>& case_weights);

  bool get_sample_with_replacement() {
    return sample_with_replacement;
  }

  const std::vector<double>& get_case_weights() {
    return case_weights;
  }

private:
  bool sample_with_replacement;
  std::vector<double> case_weights;
};

#endif //GRADIENTFOREST_SAMPLINGOPTIONS_H
