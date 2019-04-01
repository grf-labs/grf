#ifndef GRF_RCPPUTILITIES_H
#define GRF_RCPPUTILITIES_H

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "prediction/Prediction.h"
class RcppUtilities {
public:
  static Data* convert_data(Rcpp::NumericMatrix input_data,
                            Eigen::SparseMatrix<double>& sparse_input_data);

  static Rcpp::List create_prediction_object(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_prediction_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_variance_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_error_matrix(const std::vector<Prediction>& predictions);

};


#endif //GRF_RCPPUTILITIES_H
