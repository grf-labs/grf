#ifndef GRADIENTFOREST_RCPPUTILITIES_H
#define GRADIENTFOREST_RCPPUTILITIES_H

#include "forest/ForestTrainer.h"
#include "commons/globals.h"

class RcppUtilities {
public:
  static const std::string SERIALIZED_FOREST_KEY;

  static void initialize_trainer(ForestTrainer &forest_trainer,
                                 uint mtry,
                                 uint num_trees,
                                 uint num_threads,
                                 uint min_node_size,
                                 bool sample_with_replacement,
                                 double sample_fraction,
                                 std::vector<size_t> no_split_variables,
                                 uint seed,
                                 bool honesty,
                                 uint ci_group_size);

  static Rcpp::RawVector serialize_forest(const Forest& forest);

  static Forest deserialize_forest(Rcpp::RawVector input);

  static Data* convert_data(Rcpp::NumericMatrix input_data,
                            Rcpp::RawMatrix sparse_data,
                            std::vector<std::string> variable_names);

  static Rcpp::NumericMatrix create_prediction_matrix(std::vector<Prediction> predictions);

  static Rcpp::NumericMatrix create_variance_matrix(std::vector<Prediction> predictions);

};


#endif //GRADIENTFOREST_RCPPUTILITIES_H
