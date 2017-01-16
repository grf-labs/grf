#ifndef GRADIENTFOREST_FORESTTRAINERS_H
#define GRADIENTFOREST_FORESTTRAINERS_H

#include "ForestTrainer.h"

class ForestTrainers {
public:
  static ForestTrainer instrumental_trainer(Data* data,
                                            size_t outcome_index,
                                            size_t treatment_index,
                                            size_t instrument_index,
                                            double split_regularization);
  static ForestTrainer quantile_trainer(Data* data,
                                        size_t outcome_index,
                                        const std::vector<double>& quantiles);
  static ForestTrainer regression_trainer(Data* data,
                                          size_t outcome_index);
};


#endif //GRADIENTFOREST_FORESTTRAINERS_H
