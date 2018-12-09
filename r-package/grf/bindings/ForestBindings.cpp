#include <Rcpp.h>
#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestTrainer.h"

#include "commons/DefaultData.h"
#include "commons/SparseData.h"
#include "forest/ForestOptions.h"
#include "forest/Forest.h"
#include "RcppUtilities.h"
#include "serialization/ForestSerializer.h"


// [[Rcpp::export]]
Rcpp::List cpp_join_forests(const Rcpp::List forest_objects,
                       Rcpp::NumericMatrix input_data,
                       Eigen::SparseMatrix<double>& sparse_input_data) {
 Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
 std::vector<std::shared_ptr<Forest>> forest_ptrs;
 
 for (auto& forest_obj : forest_objects) {
   Forest deserialized_forest = RcppUtilities::deserialize_forest(
     static_cast<Rcpp::List>(forest_obj)[RcppUtilities::SERIALIZED_FOREST_KEY]);
   
   forest_ptrs.push_back(std::make_shared<Forest>(deserialized_forest)); 
 }

 Forest big_forest = Forest::join(forest_ptrs);
 Rcpp::List result = RcppUtilities::create_forest_object(big_forest, data);
 
 delete data;
 return result;
}
 