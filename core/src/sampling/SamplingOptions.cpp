#include "SamplingOptions.h"

SamplingOptions::SamplingOptions(size_t num_samples,
                                 bool sample_with_replacement,
                                 double sample_fraction,
                                 std::vector<double>* case_weights):
    num_samples(num_samples),
    sample_with_replacement(sample_with_replacement),
    sample_fraction(sample_fraction),
    case_weights(case_weights) {}

size_t SamplingOptions::get_num_samples() const {
  return num_samples;
}

bool SamplingOptions::get_sample_with_replacement() const {
  return sample_with_replacement;
}

double SamplingOptions::get_sample_fraction() const {
  return sample_fraction;
}

std::vector<double>* SamplingOptions::get_case_weights() const {
  return case_weights;
}
