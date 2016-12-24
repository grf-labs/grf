#include "SamplingOptions.h"

SamplingOptions::SamplingOptions(bool sample_with_replacement,
                                 double sample_fraction,
                                 std::vector<double>* case_weights,
                                 bool keep_in_bag):
    sample_with_replacement(sample_with_replacement),
    sample_fraction(sample_fraction),
    case_weights(case_weights),
    keep_in_bag(keep_in_bag) {}

bool SamplingOptions::get_sample_with_replacement() const {
  return sample_with_replacement;
}

double SamplingOptions::get_sample_fraction() const {
  return sample_fraction;
}

const std::vector<double>* SamplingOptions::get_case_weights() const {
  return case_weights;
}

bool SamplingOptions::get_keep_in_bag() const {
  return keep_in_bag;
}
