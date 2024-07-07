// PyUtilities.cpp

#include "PyUtilities.h"

py::dict PyUtilities::create_forest_object(Forest& forest, const std::vector<Prediction>& predictions) {
    py::dict forest_dict;
    forest_dict["forest"] = serialize_forest(forest);
    if (!predictions.empty()) {
        forest_dict["predictions"] = create_prediction_object(predictions);
    }
    return forest_dict;
}

py::dict PyUtilities::serialize_forest(const Forest& forest) {
    py::dict forest_dict;
    py::list trees;
    for (const auto& tree : forest.get_trees()) {
        trees.append(serialize_tree(tree));
    }
    forest_dict["trees"] = trees;
    forest_dict["num_variables"] = forest.get_num_variables();
    forest_dict["ci_group_size"] = forest.get_ci_group_size();
    return forest_dict;
}

Forest PyUtilities::deserialize_forest(const py::dict& forest_object) {
    py::list tree_list = forest_object["trees"].cast<py::list>();
    std::vector<std::unique_ptr<Tree>> trees;
    for (const auto& tree_dict : tree_list) {
        trees.push_back(deserialize_tree(tree_dict.cast<py::dict>()));
    }

    // Extract additional parameters from forest_object
    size_t num_variables = forest_object["num_variables"].cast<size_t>();
    size_t ci_group_size = forest_object["ci_group_size"].cast<size_t>();

    // Create a non-const reference to trees
    auto& trees_ref = trees;

    return Forest(trees_ref, num_variables, ci_group_size);
}

Data PyUtilities::convert_data(const py::array_t<double>& input_data) {
    py::buffer_info buf = input_data.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Number of dimensions must be 2");
    }
    size_t num_rows = buf.shape[0];
    size_t num_cols = buf.shape[1];
    double* data_ptr = static_cast<double*>(buf.ptr);
    return Data(data_ptr, num_rows, num_cols);
}

py::dict PyUtilities::create_prediction_object(const std::vector<Prediction>& predictions) {
    py::dict prediction_dict;
    prediction_dict["predictions"] = create_prediction_matrix(predictions);
    prediction_dict["variance"] = create_variance_matrix(predictions);
    prediction_dict["error"] = create_error_matrix(predictions);
    prediction_dict["excess_error"] = create_excess_error_matrix(predictions);
    return prediction_dict;
}

py::array_t<double> PyUtilities::create_prediction_matrix(const std::vector<Prediction>& predictions) {
    // Implementation depends on the structure of Prediction class
    // This is a placeholder
    return py::array_t<double>();
}

py::array_t<double> PyUtilities::create_variance_matrix(const std::vector<Prediction>& predictions) {
    // Implementation depends on the structure of Prediction class
    // This is a placeholder
    return py::array_t<double>();
}

py::array_t<double> PyUtilities::create_error_matrix(const std::vector<Prediction>& predictions) {
    // Implementation depends on the structure of Prediction class
    // This is a placeholder
    return py::array_t<double>();
}

py::array_t<double> PyUtilities::create_excess_error_matrix(const std::vector<Prediction>& predictions) {
    // Implementation depends on the structure of Prediction class
    // This is a placeholder
    return py::array_t<double>();
}

py::dict PyUtilities::serialize_tree(const std::unique_ptr<Tree>& tree) {
    py::dict tree_dict;
    tree_dict["root_node"] = tree->get_root_node();
    tree_dict["child_nodes"] = tree->get_child_nodes();
    tree_dict["leaf_samples"] = tree->get_leaf_samples();
    tree_dict["split_vars"] = tree->get_split_vars();
    tree_dict["split_values"] = tree->get_split_values();
    tree_dict["drawn_samples"] = tree->get_drawn_samples();
    tree_dict["send_missing_left"] = tree->get_send_missing_left();
    // Skipping prediction_values for now
    return tree_dict;
}

std::unique_ptr<Tree> PyUtilities::deserialize_tree(const py::dict& tree_dict) {
    size_t root_node = tree_dict["root_node"].cast<size_t>();
    std::vector<std::vector<size_t>> child_nodes = tree_dict["child_nodes"].cast<std::vector<std::vector<size_t>>>();
    std::vector<std::vector<size_t>> leaf_samples = tree_dict["leaf_samples"].cast<std::vector<std::vector<size_t>>>();
    std::vector<size_t> split_vars = tree_dict["split_vars"].cast<std::vector<size_t>>();
    std::vector<double> split_values = tree_dict["split_values"].cast<std::vector<double>>();
    std::vector<size_t> drawn_samples = tree_dict["drawn_samples"].cast<std::vector<size_t>>();
    std::vector<bool> send_missing_left = tree_dict["send_missing_left"].cast<std::vector<bool>>();
    // Using default PredictionValues for now
    PredictionValues prediction_values;
    return std::make_unique<Tree>(
        root_node, child_nodes, leaf_samples, split_vars, split_values,
        drawn_samples, send_missing_left, prediction_values
    );
}
