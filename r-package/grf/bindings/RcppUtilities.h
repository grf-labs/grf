#ifndef GRF_RCPPUTILITIES_H
#define GRF_RCPPUTILITIES_H

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestTrainer.h"

class RcppUtilities {
public:

  /**
   * Converts the provided {@link Forest} object to an R list to be returned
   * through the Rcpp bindings.
   *
   * NOTE: To converse memory, this method destructively modifies the forest
   * object by clearing out individual {@link Tree} objects. The forest cannot
   * be used once it has been passed to the method.
   */
  static Rcpp::List serialize_forest(Forest& forest);
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
