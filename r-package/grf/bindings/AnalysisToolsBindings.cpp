/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#include <Rcpp.h>
#include <vector>

#include "Eigen/Sparse"
#include "analysis/SplitFrequencyComputer.h"
#include "commons/globals.h"
#include "forest/Forest.h"
#include "prediction/collector/SampleWeightComputer.h"
#include "prediction/collector/TreeTraverser.h"

#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_split_frequencies(const Rcpp::List& forest_object,
                                              size_t max_depth) {
  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  SplitFrequencyComputer computer;
  std::vector<std::vector<size_t>> split_frequencies = computer.compute(forest, max_depth);

  size_t num_variables = forest.get_num_variables();
  Rcpp::NumericMatrix result(max_depth, num_variables);
  for (size_t depth = 0; depth < split_frequencies.size(); depth++) {
    const std::vector<size_t>& frequencies = split_frequencies.at(depth);
    for (size_t var = 0; var < num_variables; var++) {
      double frequency = frequencies[var];
        result(depth, var) = frequency;
      }
    }
  return result;
}

Eigen::SparseMatrix<double> compute_sample_weights(const Rcpp::List& forest_object,
                                                   const Rcpp::NumericMatrix& train_matrix,
                                                   const Rcpp::NumericMatrix& test_matrix,
                                                   unsigned int num_threads,
                                                   bool oob_prediction) {
  Data train_data = RcppUtilities::convert_data(train_matrix);
  Data data = RcppUtilities::convert_data(test_matrix);
  Forest forest = RcppUtilities::deserialize_forest(forest_object);
  num_threads = ForestOptions::validate_num_threads(num_threads);

  TreeTraverser tree_traverser(num_threads);
  SampleWeightComputer weight_computer;

  std::vector<std::vector<size_t>> leaf_nodes_by_tree = tree_traverser.get_leaf_nodes(forest, data, oob_prediction);
  std::vector<std::vector<bool>> trees_by_sample = tree_traverser.get_valid_trees_by_sample(forest, data, oob_prediction);

  size_t num_samples = data.get_num_rows();
  size_t num_neighbors = train_data.get_num_rows();

  // From http://eigen.tuxfamily.org/dox/group__TutorialSparse.html:
  // Filling a sparse matrix effectively
  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(num_neighbors);
  Eigen::SparseMatrix<double> result(num_samples, num_neighbors);

  for (size_t sample = 0; sample < num_samples; sample++) {
    std::unordered_map<size_t, double> weights = weight_computer.compute_weights(
        sample, forest, leaf_nodes_by_tree, trees_by_sample);
    for (auto it = weights.begin(); it != weights.end(); it++) {
      size_t neighbor = it->first;
      double weight = it->second;
      triplet_list.emplace_back(sample, neighbor, weight);
    }
  }
  result.setFromTriplets(triplet_list.begin(), triplet_list.end());

  return result;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_weights(const Rcpp::List& forest_object,
                                            const Rcpp::NumericMatrix& train_matrix,
                                            const Rcpp::NumericMatrix& test_matrix,
                                            unsigned int num_threads) {
  return compute_sample_weights(forest_object, train_matrix,
                                test_matrix, num_threads, false);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_weights_oob(const Rcpp::List& forest_object,
                                                const Rcpp::NumericMatrix& train_matrix,
                                                unsigned int num_threads) {
  return compute_sample_weights(forest_object, train_matrix,
                                train_matrix, num_threads, true);
}

// [[Rcpp::export]]
Rcpp::List merge(const Rcpp::List& forest_objects) {
 std::vector<Forest> forests;

 for (auto& forest_obj : forest_objects) {
   Forest deserialized_forest = RcppUtilities::deserialize_forest(forest_obj);
   forests.push_back(std::move(deserialized_forest));
 }

  Forest big_forest = Forest::merge(forests);
  return RcppUtilities::serialize_forest(big_forest);
}
