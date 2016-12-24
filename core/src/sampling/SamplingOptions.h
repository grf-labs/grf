#ifndef GRADIENTFOREST_SAMPLINGOPTIONS_H
#define GRADIENTFOREST_SAMPLINGOPTIONS_H


#include <string>
#include <vector>

class SamplingOptions {
public:
  SamplingOptions(bool sample_with_replacement,
                  double sample_fraction,
                  std::vector<double>* case_weights,
                  bool keep_in_bag);

  bool get_sample_with_replacement() const;
  double get_sample_fraction() const;
  std::vector<double>* get_case_weights() const;
  bool get_keep_in_bag() const;

private:
  bool sample_with_replacement;
  double sample_fraction;
  std::vector<double>* case_weights;
  bool keep_in_bag;
};

#endif //GRADIENTFOREST_SAMPLINGOPTIONS_H
