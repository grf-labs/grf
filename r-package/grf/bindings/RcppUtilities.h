#ifndef GRF_RCPPUTILITIES_H
#define GRF_RCPPUTILITIES_H

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestTrainer.h"

class RcppUtilities {
public:
  static const std::string SERIALIZED_FOREST_KEY;

  static Rcpp::List create_forest_object(const Forest& forest, Data* data);
  static Rcpp::RawVector serialize_forest(const Forest& forest);
  static Forest deserialize_forest(Rcpp::RawVector input);

  static Data* convert_data(Rcpp::NumericMatrix input_data,
                            Eigen::SparseMatrix<double>& sparse_input_data);

  static Rcpp::List create_prediction_object(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_prediction_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_variance_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_error_matrix(const std::vector<Prediction>& predictions);

};


#endif //GRF_RCPPUTILITIES_H
