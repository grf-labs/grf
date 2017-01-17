#include "SamplingOptions.h"

SamplingOptions::SamplingOptions(bool sample_with_replacement,
                                 std::vector<double>* case_weights):
    sample_with_replacement(sample_with_replacement),
    case_weights(case_weights) {}

bool SamplingOptions::get_sample_with_replacement() const {
  return sample_with_replacement;
}

std::vector<double>* SamplingOptions::get_case_weights() const {
  return case_weights;
}
