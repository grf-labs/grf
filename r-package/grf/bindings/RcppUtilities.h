#ifndef GRF_RCPPUTILITIES_H
#define GRF_RCPPUTILITIES_H

#include "forest/ForestTrainer.h"
#include "commons/globals.h"

class RcppUtilities {
public:
  static const std::string SERIALIZED_FOREST_KEY;

  static void initialize_trainer(ForestTrainer& forest_trainer,
                                 uint mtry,
                                 uint num_trees,
                                 uint num_threads,
                                 uint min_node_size,
                                 bool sample_with_replacement,
                                 double sample_fraction,
                                 const std::vector<size_t>& no_split_variables,
                                 uint seed,
                                 bool honesty,
                                 uint ci_group_size);

  static Rcpp::List create_forest_object(const Forest& forest, Data* data);
  static Rcpp::RawVector serialize_forest(const Forest& forest);
  static Forest deserialize_forest(Rcpp::RawVector input);

  static Data* convert_data(Rcpp::NumericMatrix input_data,
                            const std::vector<std::string>& variable_names);

  static Rcpp::List create_prediction_object(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_prediction_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_variance_matrix(const std::vector<Prediction>& predictions);

};


#endif //GRF_RCPPUTILITIES_H
