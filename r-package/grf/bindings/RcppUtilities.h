#ifndef GRF_RCPPUTILITIES_H
#define GRF_RCPPUTILITIES_H

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestTrainer.h"

class RcppUtilities {
public:
  static Rcpp::List serialize_forest(const Forest& forest);
  static Forest deserialize_forest(Rcpp::List forest_object);

  static Data* convert_data(Rcpp::NumericMatrix input_data,
                            Eigen::SparseMatrix<double>& sparse_input_data);

  static Rcpp::List create_prediction_object(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_prediction_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_variance_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_error_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_excess_error_matrix(const std::vector<Prediction>& predictions);

};


#endif //GRF_RCPPUTILITIES_H
