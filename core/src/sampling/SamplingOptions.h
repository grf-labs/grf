#ifndef GRADIENTFOREST_SAMPLINGOPTIONS_H
#define GRADIENTFOREST_SAMPLINGOPTIONS_H


#include <string>

class SamplingOptions {
public:
  SamplingOptions(bool sample_with_replacement,
                  double sample_fraction,
                  std::string case_weights_file);

  bool get_sample_with_replacement() const;
  double get_sample_fraction() const;
  const std::string& get_case_weights_file() const;

private:
  bool sample_with_replacement;
  double sample_fraction;
  std::string case_weights_file;
};

#endif //GRADIENTFOREST_SAMPLINGOPTIONS_H
