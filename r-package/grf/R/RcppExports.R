# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

compute_split_frequencies <- function(forest_object, max_depth) {
    .Call('_grf_compute_split_frequencies', PACKAGE = 'grf', forest_object, max_depth)
}

compute_weights <- function(forest_object, train_matrix, test_matrix, num_threads) {
    .Call('_grf_compute_weights', PACKAGE = 'grf', forest_object, train_matrix, test_matrix, num_threads)
}

compute_weights_oob <- function(forest_object, train_matrix, num_threads) {
    .Call('_grf_compute_weights_oob', PACKAGE = 'grf', forest_object, train_matrix, num_threads)
}

merge <- function(forest_objects) {
    .Call('_grf_merge', PACKAGE = 'grf', forest_objects)
}

causal_train <- function(train_matrix, outcome_index, treatment_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, reduced_form_weight, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_causal_train', PACKAGE = 'grf', train_matrix, outcome_index, treatment_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, reduced_form_weight, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

causal_predict <- function(forest_object, train_matrix, outcome_index, treatment_index, test_matrix, num_threads, estimate_variance) {
    .Call('_grf_causal_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, treatment_index, test_matrix, num_threads, estimate_variance)
}

causal_predict_oob <- function(forest_object, train_matrix, outcome_index, treatment_index, num_threads, estimate_variance) {
    .Call('_grf_causal_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, treatment_index, num_threads, estimate_variance)
}

ll_causal_predict <- function(forest_object, train_matrix, outcome_index, treatment_index, test_matrix, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance) {
    .Call('_grf_ll_causal_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, treatment_index, test_matrix, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance)
}

ll_causal_predict_oob <- function(forest_object, train_matrix, outcome_index, treatment_index, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance) {
    .Call('_grf_ll_causal_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, treatment_index, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance)
}

causal_survival_train <- function(train_matrix, causal_survival_numerator_index, causal_survival_denominator_index, treatment_index, censor_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_causal_survival_train', PACKAGE = 'grf', train_matrix, causal_survival_numerator_index, causal_survival_denominator_index, treatment_index, censor_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

causal_survival_predict <- function(forest_object, train_matrix, test_matrix, num_threads, estimate_variance) {
    .Call('_grf_causal_survival_predict', PACKAGE = 'grf', forest_object, train_matrix, test_matrix, num_threads, estimate_variance)
}

causal_survival_predict_oob <- function(forest_object, train_matrix, num_threads, estimate_variance) {
    .Call('_grf_causal_survival_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, num_threads, estimate_variance)
}

instrumental_train <- function(train_matrix, outcome_index, treatment_index, instrument_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, reduced_form_weight, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_instrumental_train', PACKAGE = 'grf', train_matrix, outcome_index, treatment_index, instrument_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, reduced_form_weight, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

instrumental_predict <- function(forest_object, train_matrix, outcome_index, treatment_index, instrument_index, test_matrix, num_threads, estimate_variance) {
    .Call('_grf_instrumental_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, treatment_index, instrument_index, test_matrix, num_threads, estimate_variance)
}

instrumental_predict_oob <- function(forest_object, train_matrix, outcome_index, treatment_index, instrument_index, num_threads, estimate_variance) {
    .Call('_grf_instrumental_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, treatment_index, instrument_index, num_threads, estimate_variance)
}

multi_causal_train <- function(train_matrix, outcome_index, treatment_index, sample_weight_index, use_sample_weights, gradient_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_multi_causal_train', PACKAGE = 'grf', train_matrix, outcome_index, treatment_index, sample_weight_index, use_sample_weights, gradient_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, stabilize_splits, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

multi_causal_predict <- function(forest_object, train_matrix, test_matrix, num_outcomes, num_treatments, num_threads, estimate_variance) {
    .Call('_grf_multi_causal_predict', PACKAGE = 'grf', forest_object, train_matrix, test_matrix, num_outcomes, num_treatments, num_threads, estimate_variance)
}

multi_causal_predict_oob <- function(forest_object, train_matrix, num_outcomes, num_treatments, num_threads, estimate_variance) {
    .Call('_grf_multi_causal_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, num_outcomes, num_treatments, num_threads, estimate_variance)
}

multi_regression_train <- function(train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_multi_regression_train', PACKAGE = 'grf', train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

multi_regression_predict <- function(forest_object, train_matrix, test_matrix, num_outcomes, num_threads) {
    .Call('_grf_multi_regression_predict', PACKAGE = 'grf', forest_object, train_matrix, test_matrix, num_outcomes, num_threads)
}

multi_regression_predict_oob <- function(forest_object, train_matrix, num_outcomes, num_threads) {
    .Call('_grf_multi_regression_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, num_outcomes, num_threads)
}

probability_train <- function(train_matrix, outcome_index, sample_weight_index, use_sample_weights, num_classes, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_probability_train', PACKAGE = 'grf', train_matrix, outcome_index, sample_weight_index, use_sample_weights, num_classes, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

probability_predict <- function(forest_object, train_matrix, outcome_index, num_classes, test_matrix, num_threads, estimate_variance) {
    .Call('_grf_probability_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, num_classes, test_matrix, num_threads, estimate_variance)
}

probability_predict_oob <- function(forest_object, train_matrix, outcome_index, num_classes, num_threads, estimate_variance) {
    .Call('_grf_probability_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, num_classes, num_threads, estimate_variance)
}

quantile_train <- function(quantiles, regression_splitting, train_matrix, outcome_index, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_quantile_train', PACKAGE = 'grf', quantiles, regression_splitting, train_matrix, outcome_index, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

quantile_predict <- function(forest_object, quantiles, train_matrix, outcome_index, test_matrix, num_threads) {
    .Call('_grf_quantile_predict', PACKAGE = 'grf', forest_object, quantiles, train_matrix, outcome_index, test_matrix, num_threads)
}

quantile_predict_oob <- function(forest_object, quantiles, train_matrix, outcome_index, num_threads) {
    .Call('_grf_quantile_predict_oob', PACKAGE = 'grf', forest_object, quantiles, train_matrix, outcome_index, num_threads)
}

regression_train <- function(train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed) {
    .Call('_grf_regression_train', PACKAGE = 'grf', train_matrix, outcome_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, compute_oob_predictions, num_threads, seed, legacy_seed)
}

regression_predict <- function(forest_object, train_matrix, outcome_index, test_matrix, num_threads, estimate_variance) {
    .Call('_grf_regression_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, test_matrix, num_threads, estimate_variance)
}

regression_predict_oob <- function(forest_object, train_matrix, outcome_index, num_threads, estimate_variance) {
    .Call('_grf_regression_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, num_threads, estimate_variance)
}

ll_regression_train <- function(train_matrix, outcome_index, ll_split_lambda, ll_split_weight_penalty, ll_split_variables, ll_split_cutoff, overall_beta, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, num_threads, seed, legacy_seed) {
    .Call('_grf_ll_regression_train', PACKAGE = 'grf', train_matrix, outcome_index, ll_split_lambda, ll_split_weight_penalty, ll_split_variables, ll_split_cutoff, overall_beta, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, ci_group_size, alpha, imbalance_penalty, clusters, samples_per_cluster, num_threads, seed, legacy_seed)
}

ll_regression_predict <- function(forest_object, train_matrix, outcome_index, test_matrix, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance) {
    .Call('_grf_ll_regression_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, test_matrix, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance)
}

ll_regression_predict_oob <- function(forest_object, train_matrix, outcome_index, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance) {
    .Call('_grf_ll_regression_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, ll_lambda, ll_weight_penalty, linear_correction_variables, num_threads, estimate_variance)
}

survival_train <- function(train_matrix, outcome_index, censor_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, alpha, num_failures, clusters, samples_per_cluster, compute_oob_predictions, prediction_type, num_threads, seed, legacy_seed) {
    .Call('_grf_survival_train', PACKAGE = 'grf', train_matrix, outcome_index, censor_index, sample_weight_index, use_sample_weights, mtry, num_trees, min_node_size, sample_fraction, honesty, honesty_fraction, honesty_prune_leaves, alpha, num_failures, clusters, samples_per_cluster, compute_oob_predictions, prediction_type, num_threads, seed, legacy_seed)
}

survival_predict <- function(forest_object, train_matrix, outcome_index, censor_index, sample_weight_index, use_sample_weights, prediction_type, test_matrix, num_threads, num_failures) {
    .Call('_grf_survival_predict', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, censor_index, sample_weight_index, use_sample_weights, prediction_type, test_matrix, num_threads, num_failures)
}

survival_predict_oob <- function(forest_object, train_matrix, outcome_index, censor_index, sample_weight_index, use_sample_weights, prediction_type, num_threads, num_failures) {
    .Call('_grf_survival_predict_oob', PACKAGE = 'grf', forest_object, train_matrix, outcome_index, censor_index, sample_weight_index, use_sample_weights, prediction_type, num_threads, num_failures)
}

