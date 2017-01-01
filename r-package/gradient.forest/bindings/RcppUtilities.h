#ifndef GRADIENTFOREST_RCPPUTILITIES_H
#define GRADIENTFOREST_RCPPUTILITIES_H

#include "ForestModel.h"
#include "globals.h"

class RcppUtilities {
public:
  static const std::string SERIALIZED_FOREST_KEY;

  static void initializeForestModel(ForestModel *forest_model,
                                    uint mtry,
                                    uint num_trees,
                                    uint num_threads,
                                    uint min_node_size,
                                    bool sample_with_replacement,
                                    double sample_fraction);

  static Rcpp::RawVector serialize_forest(Forest *forest);

  static Forest* deserialize_forest(Rcpp::RawVector input);

  static Data* convert_data(Rcpp::NumericMatrix input_data,
                            Rcpp::RawMatrix sparse_data,
                            std::vector<std::string> variable_names);

  static Rcpp::NumericMatrix create_prediction_matrix(std::vector<std::vector<double>> predictions,
                                                      size_t prediction_length);
};


#endif //GRADIENTFOREST_RCPPUTILITIES_H
