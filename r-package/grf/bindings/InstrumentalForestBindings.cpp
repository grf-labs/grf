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
                              std::vector<std::string> variable_names,
                              unsigned int mtry,
                              unsigned int num_trees,
                              bool verbose,
                              unsigned int num_threads,
                              unsigned int min_node_size,
                              bool sample_with_replacement,
                              bool keep_inbag,
                              double sample_fraction,
                              unsigned int seed,
                              bool honesty,
                              unsigned int ci_group_size,
                              double split_regularization,
                              double alpha,
                              double lambda,
                              bool downweight_penalty) {
  ForestTrainer trainer = lambda > 0
    ? ForestTrainers::regularized_instrumental_trainer(outcome_index - 1,
                                                       treatment_index - 1,
                                                       instrument_index - 1,
                                                       split_regularization,
                                                       lambda,
                                                       downweight_penalty)
    : ForestTrainers::instrumental_trainer(outcome_index - 1,
                                           treatment_index - 1,
                                           instrument_index - 1,
                                           split_regularization,
                                           alpha);

  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);
  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size,
                        honesty, sample_with_replacement, num_threads, seed);

  Forest forest = trainer.train(data, options);

  Rcpp::List result = RcppUtilities::create_forest_object(forest, data);
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List instrumental_predict(Rcpp::List forest_object,
                                Rcpp::NumericMatrix input_data,
                                Eigen::SparseMatrix<double> sparse_input_data,
                                std::vector <std::string> variable_names,
                                unsigned int num_threads,
                                unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);
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
                                    std::vector <std::string> variable_names,
                                    unsigned int num_threads,
                                    unsigned int ci_group_size) {
  Data* data = RcppUtilities::convert_data(input_data, sparse_input_data, variable_names);
  Forest forest = RcppUtilities::deserialize_forest(
      forest_object[RcppUtilities::SERIALIZED_FOREST_KEY]);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads, ci_group_size);
  std::vector<Prediction> predictions = predictor.predict_oob(forest, data);

  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);
  delete data;
  return result;
}
