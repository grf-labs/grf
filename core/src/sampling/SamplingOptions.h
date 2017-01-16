#ifndef GRADIENTFOREST_SAMPLINGOPTIONS_H
#define GRADIENTFOREST_SAMPLINGOPTIONS_H


#include <string>
#include <vector>

class SamplingOptions {
public:
  SamplingOptions(size_t num_samples,
                  bool sample_with_replacement,
                  double sample_fraction,
                  std::vector<double>* case_weights);

  size_t get_num_samples() const;
  bool get_sample_with_replacement() const;
  double get_sample_fraction() const;
  std::vector<double>* get_case_weights() const;

private:
  size_t num_samples;
  bool sample_with_replacement;
  double sample_fraction;
  std::vector<double>* case_weights;
};

#endif //GRADIENTFOREST_SAMPLINGOPTIONS_H
