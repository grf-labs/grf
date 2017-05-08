/*-------------------------------------------------------------------------------
  This file is part of gradient-forest.

  gradient-forest is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gradient-forest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with gradient-forest. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List custom_train(Rcpp::NumericMatrix input_data,
                        size_t outcome_index,
                        Rcpp::RawMatrix sparse_data,
                        std::vector <std::string> variable_names,
                        unsigned int mtry,
                        unsigned int num_trees,
                        bool verbose,
                        unsigned int num_threads,
                        unsigned int min_node_size,
                        bool sample_with_replacement,
                        bool keep_inbag,
                        double sample_fraction,
                        std::vector<size_t> no_split_variables,
                        unsigned int seed,
                        bool honesty,
                        unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  ForestTrainer trainer = ForestTrainers::custom_trainer(data, outcome_index);
  RcppUtilities::initialize_trainer(trainer, mtry, num_trees, num_threads, min_node_size,
                                    sample_with_replacement, sample_fraction, no_split_variables, seed, honesty, ci_group_size);
  Forest forest = trainer.train(data);

  Rcpp::List result;
  Rcpp::RawVector serialized_forest = RcppUtilities::serialize_forest(forest);
  result.push_back(serialized_forest, RcppUtilities::SERIALIZED_FOREST_KEY);
  result.push_back(forest.get_trees().size(), "num.trees");

  delete data;

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix custom_predict(Rcpp::List forest,
                                   Rcpp::NumericMatrix input_data,
                                   Rcpp::RawMatrix sparse_data,
                                   std::vector<std::string> variable_names,
                                   unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::custom_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix custom_predict_oob(Rcpp::List forest,
                                       Rcpp::NumericMatrix input_data,
                                       Rcpp::RawMatrix sparse_data,
                                       std::vector<std::string> variable_names,
                                       unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest deserialized_forest = RcppUtilities::deserialize_forest(
      forest[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::custom_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(deserialized_forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}
