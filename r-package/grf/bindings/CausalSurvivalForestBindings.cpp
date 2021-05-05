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
#include <vector>

#include "commons/globals.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

using namespace grf;

// [[Rcpp::export]]
Rcpp::List causal_survival_train(const Rcpp::NumericMatrix& train_matrix,
                                 size_t causal_survival_numerator_index,
                                 size_t causal_survival_denominator_index,
                                 size_t treatment_index,
                                 size_t censor_index,
                                 size_t sample_weight_index,
                                 bool use_sample_weights,
                                 unsigned int mtry,
                                 unsigned int num_trees,
                                 unsigned int min_node_size,
                                 double sample_fraction,
                                 bool honesty,
                                 double honesty_fraction,
                                 bool honesty_prune_leaves,
                                 size_t ci_group_size,
                                 double alpha,
                                 double imbalance_penalty,
                                 bool stabilize_splits,
                                 const std::vector<size_t>& clusters,
                                 unsigned int samples_per_cluster,
                                 bool compute_oob_predictions,
                                 unsigned int num_threads,
                                 unsigned int seed) {
  ForestTrainer trainer = causal_survival_trainer(stabilize_splits);

  Data data = RcppUtilities::convert_data(train_matrix);
  data.set_causal_survival_numerator_index(causal_survival_numerator_index);
  data.set_causal_survival_denominator_index(causal_survival_denominator_index);
  data.set_treatment_index(treatment_index);
  data.set_instrument_index(treatment_index);
  data.set_censor_index(censor_index);
  if (use_sample_weights) {
    data.set_weight_index(sample_weight_index);
  }

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
  Forest forest = trainer.train(data, options);

  std::vector<Prediction> predictions;
  if (compute_oob_predictions) {
    ForestPredictor predictor = causal_survival_predictor(num_threads);
    predictions = predictor.predict_oob(forest, data, false);
  }

  return RcppUtilities::create_forest_object(forest, predictions);
}

// [[Rcpp::export]]
Rcpp::List causal_survival_predict(const Rcpp::List& forest_object,
                                   const Rcpp::NumericMatrix& train_matrix,
                                   const Rcpp::NumericMatrix& test_matrix,
                                   unsigned int num_threads,
                                   bool estimate_variance) {
  Data train_data = RcppUtilities::convert_data(train_matrix);
  Data data = RcppUtilities::convert_data(test_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = causal_survival_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(forest, train_data, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}

// [[Rcpp::export]]
Rcpp::List causal_survival_predict_oob(const Rcpp::List& forest_object,
                                       const Rcpp::NumericMatrix& train_matrix,
                                       unsigned int num_threads,
                                       bool estimate_variance) {
  Data data = RcppUtilities::convert_data(train_matrix);

  Forest forest = RcppUtilities::deserialize_forest(forest_object);

  ForestPredictor predictor = causal_survival_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  return result;
}
