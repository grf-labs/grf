#' Quantile forest
#' 
#' Trains a regression forest that can be used to estimate
#' quantiles of the conditional distribution of Y given X = x.
#'
#' @param X The covariates used in the quantile regression.
#' @param Y The outcome.
#' @param quantiles Vector of quantiles used to calibrate the forest.
#' @param regression.splitting Whether to use regression splits when growing trees instead
#'                             of specialized splits based on the quantiles (the default).
#'                             Setting this flag to true corresponds to the approach to
#'                             quantile forests from Meinshausen (2006).
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty is used, these subsamples will
#'                        further be cut in half.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#' @param seed The seed for the C++ random number generator.
#' @param honesty Whether or not honest splitting (i.e., sub-sample splitting) should be used.
#' @param alpha Maximum imbalance of a split.
#'
#' @return A trained quantile forest object.
#'
#' @examples \dontrun{
#' # Generate data.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' Y = X[,1] * rnorm(n)
#' 
#' # Train a quantile forest.
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#' 
#' # Make predictions.
#' q.hat = predict(q.forest, X.test)
#' 
#' # Make predictions for different quantiles than those used in training.
#' q.hat = predict(q.forest, X.test, quantiles=c(0.1, 0.9))
#' 
#' # Train a quantile forest using regression splitting instead of quantile-based
#' # splits, emulating the approach in Meinshausen (2006).
#' meins.forest = quantile_forest(X, Y, regression.splitting=TRUE)
#' 
#' # Make predictions for the desired quantiles.
#' q.hat = predict(meins.forest, X.test, quantiles=c(0.1, 0.5, 0.9))
#' }
#' 
#' @export
quantile_forest <- function(X, Y, quantiles = c(0.1, 0.5, 0.9), regression.splitting = FALSE,
                            sample.fraction = 0.5, mtry = NULL, num.trees = 2000,
                            num.threads = NULL, min.node.size = NULL, seed = NULL, alpha = 0.05,
                            honesty = TRUE) {
    
    if (!is.numeric(quantiles) | length(quantiles) < 1) {
        stop("Error: Must provide numeric quantiles")
    } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
        stop("Error: Quantiles must be in (0, 1)")
    }
    
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    
    mtry <- validate_mtry(mtry, X)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    
    no.split.variables <- numeric(0)
    sample.with.replacement <- FALSE
    verbose <- FALSE
    keep.inbag <- FALSE
    
    data <- create_data_matrices(X, Y)
    variable.names <- c(colnames(X), "outcome")
    outcome.index <- ncol(X) + 1

    ci.group.size <- 1
    
    forest <- quantile_train(quantiles, regression.splitting, data$default, data$sparse, outcome.index,
        variable.names, mtry, num.trees, verbose, num.threads, min.node.size, sample.with.replacement,
        keep.inbag, sample.fraction, no.split.variables, seed, honesty, ci.group.size, alpha)
    
    forest[["X.orig"]] <- X
    class(forest) <- c("quantile_forest", "grf")
    forest
}

#' Predict with a quantile forest
#' 
#' Gets estimates of the conditional quantiles of Y given X using a trained forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param quantiles Vector of quantiles at which estimates are required.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Predictions at each test point for each desired quantile.
#'
#' @examples \dontrun{
#' # Train a quantile forest.
#' n = 50; p = 10
#' X = matrix(rnorm(n*p), n, p)
#' Y = X[,1] * rnorm(n)
#' q.forest = quantile_forest(X, Y, quantiles=c(0.1, 0.5, 0.9))
#'
#' # Predict on out-of-bag training samples.
#' q.pred = predict(q.forest)
#'
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[,1] = seq(-2, 2, length.out = 101)
#' q.pred = predict(q.forest, X.test)
#' }
#'
#' @export
predict.quantile_forest <- function(object,
                                    newdata = NULL,
                                    quantiles = c(0.1, 0.5, 0.9),
                                    num.threads = NULL, ...) {
    
    if (!is.numeric(quantiles) | length(quantiles) < 1) {
        stop("Error: Must provide numeric quantiles")
    } else if (min(quantiles) <= 0 | max(quantiles) >= 1) {
        stop("Error: Quantiles must be in (0, 1)")
    }
    
    num.threads <- validate_num_threads(num.threads)
    variable.names <- character(0)
    
    forest.short <- object[-which(names(object) == "X.orig")]
    
    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata, NA)
        quantile_predict(forest.short, quantiles, data$default, data$sparse,
                         variable.names, num.threads)
    } else {
        data <- create_data_matrices(object[["X.orig"]], NA)
        quantile_predict_oob(forest.short, quantiles, data$default, data$sparse,
                             variable.names, num.threads)
    }
}
