#include <map>
#include <Rcpp.h>
#include <sstream>
#include <vector>

#include "commons/globals.h"
#include "Eigen/Sparse"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainers.h"
#include "RcppUtilities.h"

// [[Rcpp::export]]
Rcpp::List instrumental_train(Rcpp::NumericMatrix input_data,
                              Eigen::SparseMatrix<double> sparse_input_data,
                              size_t outcome_index,
                              size_t treatment_index,
                              size_t instrument_index,
                              unsigned int mtry,
                              unsigned int num_trees,
                              unsigned int num_threads,
                              unsigned int min_node_size,
                              double sample_fraction,
                              unsigned int seed,
                              bool honesty,
                              double honesty_fraction,
                              unsigned int ci_group_size,
                              double reduced_form_weight,
                              double alpha,
                              bool imbalance_penalty,
                              bool stabilize_splits,
                              std::vector<size_t> clusters,
                              unsigned int samples_per_cluster) {
  ForestTrainer trainer = ForestTrainers::instrumental_trainer(
      outcome_index - 1,
      treatment_index - 1,
      instrument_index - 1,
      reduced_form_weight,
      stabilize_splits);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List instrumental_predict(Rcpp::List forest_object,
                                Rcpp::NumericMatrix input_data,
                                Eigen::SparseMatrix<double> sparse_input_data,
                                unsigned int num_threads,
                                unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List instrumental_predict_oob(Rcpp::List forest_object,
                                    Rcpp::NumericMatrix input_data,
                                    Eigen::SparseMatrix<double> sparse_input_data,
                                    unsigned int num_threads,
                                    unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}
