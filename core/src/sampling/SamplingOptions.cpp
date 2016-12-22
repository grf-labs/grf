#include "SamplingOptions.h"

SamplingOptions::SamplingOptions(bool sample_with_replacement,
                                 double sample_fraction,
                                 std::string case_weights_file):
    sample_with_replacement(sample_with_replacement),
    sample_fraction(sample_fraction),
    case_weights_file(case_weights_file) {}

bool SamplingOptions::get_sample_with_replacement() const {
  return sample_with_replacement;
}

double SamplingOptions::get_sample_fraction() const {
  return sample_fraction;
}

const std::string &SamplingOptions::get_case_weights_file() const {
  return case_weights_file;
}
