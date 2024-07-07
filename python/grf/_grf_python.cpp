// _grf_python.cpp

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "forest/ForestTrainers.h"
#include "forest/ForestPredictors.h"
#include "PyUtilities.h"

namespace py = pybind11;
using namespace grf;

py::object regression_train(
    py::array_t<double> train_matrix,
    size_t outcome_index,
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
    std::vector<size_t> clusters,
    unsigned int samples_per_cluster,
    bool compute_oob_predictions,
    unsigned int num_threads,
    unsigned int seed)
{
    ForestTrainer trainer = regression_trainer();

    Data data = PyUtilities::convert_data(train_matrix);
    data.set_outcome_index(outcome_index);
    if (use_sample_weights) {
        data.set_weight_index(sample_weight_index);
    }

    ForestOptions options(num_trees, ci_group_size, sample_fraction, mtry, min_node_size, honesty,
        honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, num_threads, seed, clusters, samples_per_cluster);
    Forest forest = trainer.train(data, options);

    std::vector<Prediction> predictions;
    if (compute_oob_predictions) {
        ForestPredictor predictor = regression_predictor(num_threads);
        predictions = predictor.predict_oob(forest, data, false);
    }

    return PyUtilities::create_forest_object(forest, predictions);
}

py::dict regression_predict(
    py::object forest_object,
    py::array_t<double> train_matrix,
    size_t outcome_index,
    py::array_t<double> test_matrix,
    unsigned int num_threads,
    bool estimate_variance)
{
    Data train_data = PyUtilities::convert_data(train_matrix);
    train_data.set_outcome_index(outcome_index);

    Data data = PyUtilities::convert_data(test_matrix);
    Forest forest = PyUtilities::deserialize_forest(forest_object);

    // Add debug information
    std::cout << "Number of trees in forest: " << forest.get_trees().size() << std::endl;
    std::cout << "Number of variables: " << forest.get_num_variables() << std::endl;
    std::cout << "Train data dimensions: " << train_data.get_num_rows() << "x" << train_data.get_num_cols() << std::endl;
    std::cout << "Test data dimensions: " << data.get_num_rows() << "x" << data.get_num_cols() << std::endl;

    ForestPredictor predictor = regression_predictor(num_threads);
    std::vector<Prediction> predictions;
    try {
        predictions = predictor.predict(forest, train_data, data, estimate_variance);
    } catch (const std::exception& e) {
        std::cerr << "Error in predict: " << e.what() << std::endl;
        throw;
    }

    return PyUtilities::create_prediction_object(predictions);
}

py::dict regression_predict_oob(
    py::object forest_object,
    py::array_t<double> train_matrix,
    size_t outcome_index,
    unsigned int num_threads,
    bool estimate_variance)
{
    Data data = PyUtilities::convert_data(train_matrix);
    data.set_outcome_index(outcome_index);

    Forest forest = PyUtilities::deserialize_forest(forest_object);

    ForestPredictor predictor = regression_predictor(num_threads);
    std::vector<Prediction> predictions = predictor.predict_oob(forest, data, estimate_variance);

    return PyUtilities::create_prediction_object(predictions);
}

PYBIND11_MODULE(_grf_python, m) {
    m.def("regression_train", &regression_train, "Train a regression forest");
    m.def("regression_predict", &regression_predict, "Predict using a regression forest");
    m.def("regression_predict_oob", &regression_predict_oob, "Predict OOB using a regression forest");
}
