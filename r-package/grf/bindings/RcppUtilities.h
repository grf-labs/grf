/*-------------------------------------------------------------------------------
  Copyright (c) 2024 GRF Contributors.

  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_RCPPUTILITIES_H
#define GRF_RCPPUTILITIES_H

#include "commons/globals.h"
#include "forest/ForestTrainer.h"

using namespace grf;

class RcppUtilities {
public:

  /**
   * Converts the provided {@link Forest} object and OOB predictions to an R list
   * to be returned through the Rcpp bindings. The provided predictions vector can
   * be present if OOB predictions were not requested as part of training.
   *
   * NOTE: To conserve memory, this method destructively modifies the forest
   * object by clearing out individual {@link Tree} objects. The forest cannot
   * be used once it has been passed to the method.
   */
  static Rcpp::List create_forest_object(Forest& forest,
                                         const std::vector<Prediction>& predictions);

  static Rcpp::List serialize_forest(Forest& forest);
  static Forest deserialize_forest(const Rcpp::List& forest_object);

  static Data convert_data(const Rcpp::NumericMatrix& input_data);

  static Rcpp::List create_prediction_object(const std::vector<Prediction>& predictions);
  static void add_predictions(Rcpp::List& output,
                              const std::vector<Prediction>& predictions);

  static Rcpp::NumericMatrix create_prediction_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_variance_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_error_matrix(const std::vector<Prediction>& predictions);
  static Rcpp::NumericMatrix create_excess_error_matrix(const std::vector<Prediction>& predictions);

};


#endif //GRF_RCPPUTILITIES_H
