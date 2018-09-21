#' Custom forest
#' 
#' Trains a custom forest model.
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting).
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds 
#'                         to set J1 in the notation of the paper. When using the defaults (honesty = TRUE and 
#'                         honesty.fraction = NULL), half of the data will be used for determining splits
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param seed The seed for the C++ random number generator.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples_per_cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster. Must be less than the size of the smallest cluster. If set to NULL
#'                            software will set this value to the size of the smallest cluster.
#'
#' @return A trained regression forest object.
#'
#' @examples \dontrun{
#' # Train a custom forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' c.forest = custom_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#' }
#'
#' @export
custom_forest <- function(X, Y, sample.fraction = 0.5, mtry = NULL, 
    num.trees = 2000, num.threads = NULL, min.node.size = NULL, honesty = TRUE,
    honesty.fraction = NULL, alpha = 0.05, imbalance.penalty = 0.0, seed = NULL,
    clusters = NULL, samples_per_cluster = NULL) {

    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    
    mtry <- validate_mtry(mtry, X)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)
    honesty.fraction <- validate_honesty_fraction(honesty.fraction, honesty)
    
    no.split.variables <- numeric(0)
    
    data <- create_data_matrices(X, Y)
    outcome.index <- ncol(X) + 1
    ci.group.size <- 1
    
    forest <- custom_train(data$default, data$sparse, outcome.index, mtry,num.trees, num.threads,
        min.node.size, sample.fraction, seed, honesty, coerce_honesty_fraction(honesty.fraction),
        ci.group.size, alpha, imbalance.penalty, clusters, samples_per_cluster)
    
    forest[["X.orig"]] <- X
    class(forest) <- c("custom_forest", "grf")
    forest
}

#' Predict with a custom forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions.
#'
#' @examples \dontrun{
#' # Train a custom forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' c.forest = custom_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' c.pred = predict(c.forest, X.test)
#' }
#'
#' @method predict custom_forest
#' @export
predict.custom_forest <- function(object, newdata = NULL, num.threads = NULL, ...) {
    
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
        
    forest.short <- object[-which(names(object) == "X.orig")]
    
    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata)
        custom_predict(forest.short, data$default, data$sparse, num.threads)
    } else {
        data <- create_data_matrices(object[["X.orig"]])
        custom_predict_oob(forest.short, data$default, data$sparse, num.threads)
    }
}
