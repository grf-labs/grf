#include "SamplingOptions.h"

SamplingOptions::SamplingOptions(bool sample_with_replacement,
                                 const std::vector<double>& case_weights):
    sample_with_replacement(sample_with_replacement),
    case_weights(case_weights) {}
