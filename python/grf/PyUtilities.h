// PyUtilities.h

#ifndef PY_UTILITIES_H
#define PY_UTILITIES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <memory>
#include <vector>
#include "forest/Forest.h"
#include "tree/Tree.h"
#include "commons/Data.h"
#include "prediction/Prediction.h"

namespace py = pybind11;
using namespace grf;

class PyUtilities {
public:
    static py::dict create_forest_object(Forest& forest, const std::vector<Prediction>& predictions);
    static py::dict serialize_forest(const Forest& forest);
    static Forest deserialize_forest(const py::dict& forest_object);
    static Data convert_data(const py::array_t<double>& input_data);
    static py::dict create_prediction_object(const std::vector<Prediction>& predictions);
    static py::array_t<double> create_prediction_matrix(const std::vector<Prediction>& predictions);
    static py::array_t<double> create_variance_matrix(const std::vector<Prediction>& predictions);
    static py::array_t<double> create_error_matrix(const std::vector<Prediction>& predictions);
    static py::array_t<double> create_excess_error_matrix(const std::vector<Prediction>& predictions);

private:
    static py::dict serialize_tree(const std::unique_ptr<Tree>& tree);
    static std::unique_ptr<Tree> deserialize_tree(const py::dict& tree_dict);
};

#endif // PY_UTILITIES_H
