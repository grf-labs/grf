/*-------------------------------------------------------------------------------
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
                        unsigned int ci_group_size,
                        double alpha) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);

  ForestTrainer trainer = ForestTrainers::custom_trainer(data,
          outcome_index - 1,
          alpha);
  RcppUtilities::initialize_trainer(trainer, mtry, num_trees, num_threads, min_node_size,
                                    sample_with_replacement, sample_fraction, no_split_variables, seed, honesty, ci_group_size);
  Forest forest = trainer.train(data);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix custom_predict(Rcpp::List forest_object,
                                   Rcpp::NumericMatrix input_data,
                                   Rcpp::RawMatrix sparse_data,
                                   std::vector<std::string> variable_names,
                                   unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::custom_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix custom_predict_oob(Rcpp::List forest_object,
                                       Rcpp::NumericMatrix input_data,
                                       Rcpp::RawMatrix sparse_data,
                                       std::vector<std::string> variable_names,
                                       unsigned int num_threads) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_data, variable_names);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::custom_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);
  Rcpp::NumericMatrix result = RcppUtilities::create_prediction_matrix(predictions);

  delete data;
  return result;
}
