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
SEXP instrumental_train(Rcpp::NumericMatrix train_matrix,
                              Eigen::SparseMatrix<double> sparse_train_matrix,
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
                              size_t ci_group_size,
                              double reduced_form_weight,
                              double alpha,
                              double imbalance_penalty,
                              bool stabilize_splits,
                              std::vector<size_t> clusters,
                              unsigned int samples_per_cluster) {
  ForestTrainer trainer = ForestTrainers::instrumental_trainer(reduced_form_weight, stabilize_splits);

  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_treatment_index(treatment_index - 1);
  data->set_instrument_index(instrument_index - 1);
  data->sort();

  ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
      honesty_fraction, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);

  Forest forest = trainer.train(data, options);

  delete data;
  return Rcpp::XPtr<Forest>(new Forest(forest));
}

// [[Rcpp::export]]
Rcpp::List instrumental_predict(SEXP xp,
                                Rcpp::NumericMatrix train_matrix,
                                Eigen::SparseMatrix<double> sparse_train_matrix,
                                size_t outcome_index,
                                size_t treatment_index,
                                size_t instrument_index,
                                Rcpp::NumericMatrix test_matrix,
                                Eigen::SparseMatrix<double> sparse_test_matrix,
                                unsigned int num_threads,
                                bool estimate_variance) {
  Rcpp::XPtr<Forest> forest(xp);
  Data* train_data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  train_data->set_outcome_index(outcome_index - 1);
  train_data->set_treatment_index(treatment_index - 1);
  train_data->set_instrument_index(instrument_index - 1);
  Data* data = RcppUtilities::convert_data(test_matrix, sparse_test_matrix);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict(*forest, train_data, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete train_data;
  delete data;
  return result;
}

// [[Rcpp::export]]
Rcpp::List instrumental_predict_oob(SEXP xp,
                                    Rcpp::NumericMatrix train_matrix,
                                    Eigen::SparseMatrix<double> sparse_train_matrix,
                                    size_t outcome_index,
                                    size_t treatment_index,
                                    size_t instrument_index,
                                    unsigned int num_threads,
                                    bool estimate_variance) {
  Rcpp::XPtr<Forest> forest(xp);
  Data* data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);
  data->set_outcome_index(outcome_index - 1);
  data->set_treatment_index(treatment_index - 1);
  data->set_instrument_index(instrument_index - 1);

  ForestPredictor predictor = ForestPredictors::instrumental_predictor(num_threads);
  std::vector<Prediction> predictions = predictor.predict_oob(*forest, data, estimate_variance);
  Rcpp::List result = RcppUtilities::create_prediction_object(predictions);

  delete data;
  return result;
}
