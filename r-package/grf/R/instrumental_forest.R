#' Intrumental forest
#' 
#' Trains an instrumental forest that can be used to estimate
#' conditional local average treatment effects tau(X) identified
#' using instruments. Formally, the forest estimates
#' tau(X) = Cov[Y, Z | X = x] / Cov[W, Z | X = x].
#' Note that when the instrument Z and treatment assignment W
#' coincide, an instrumental forest is equivalent to a causal forest.
#'
#' @param X The covariates used in the instrumental regression.
#' @param Y The outcome.
#' @param W The treatment assignment (may be binary or real).
#' @param Z The instrument (may be binary or real).
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty is used, these subsamples will
#'                        further be cut in half.
#' @param mtry Number of variables tried for each split.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param min.node.size Minimum number of observations in each tree leaf.
#' @param honesty Should honest splitting (i.e., sub-sample splitting) be used?
#' @param ci.group.size The forst will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param precompute.nuisance Should we first run regression forests to estimate
#'                            y(x) = E[Y|X=x], w(x) = E[W|X=x] and z(x) = E[Z|X=x],
#'                            and then run an instrumental forest on the residuals?
#'                            This approach is recommended, computational resources
#'                            permitting.
#' @param split.regularization Whether splits should be regularized towards a naive
#'                             splitting criterion that ignores the instrument (and
#'                             instead emulates a causal forest).
#' @param alpha Maximum imbalance of a split.
#' @param seed The seed of the c++ random number generator.
#'
#' @return A trained instrumental forest object.
#' @export
instrumental_forest <- function(X, Y, W, Z, sample.fraction = 0.5, mtry = ceiling(2*ncol(X)/3), 
                                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                honesty = TRUE, ci.group.size = 2, precompute.nuisance = TRUE,
                                split.regularization = 0, alpha = 0.05, seed = NULL) {
    
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    if(length(W) != nrow(X)) { stop("W has incorrect length.") }
    if(length(Z) != nrow(X)) { stop("Z has incorrect length.") }
    
    mtry <- validate_mtry(mtry)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    
    sparse.data <- as.matrix(0)
    no.split.variables <- numeric(0)
    sample.with.replacement <- FALSE
    verbose <- FALSE
    keep.inbag <- FALSE
    
    if (!is.numeric(split.regularization) | split.regularization < 0 | split.regularization > 1) {
        stop("Error: Invalid value for split.regularization. Please give a value in [0,1].")
    }
    
    if (!precompute.nuisance) {
        
        input.data <- as.matrix(cbind(X, Y, W, Z))
        
    } else {
        
        forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction, mtry = mtry, 
                                      num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                      honesty = TRUE, seed = seed, ci.group.size = 1)
        Y.hat = predict(forest.Y)$predictions
        
        forest.W <- regression_forest(X, W, sample.fraction = sample.fraction, mtry = mtry, 
                                      num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                      honesty = TRUE, seed = seed, ci.group.size = 1)
        W.hat = predict(forest.W)$predictions
        
        forest.Z <- regression_forest(X, Z, sample.fraction = sample.fraction, mtry = mtry, 
                                      num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, 
                                      honesty = TRUE, seed = seed, ci.group.size = 1)
        Z.hat = predict(forest.Z)$predictions
        
        input.data <- as.matrix(cbind(X, Y - Y.hat, W - W.hat, Z - Z.hat))
        
    }
    
    variable.names <- c(colnames(X), "outcome", "treatment", "instrument")
    outcome.index <- ncol(X) + 1
    treatment.index <- ncol(X) + 2
    instrument.index <- ncol(X) + 3
    
    forest <- instrumental_train(input.data, outcome.index, treatment.index,
        instrument.index, sparse.data, variable.names, mtry, num.trees, verbose,
        num.threads, min.node.size, sample.with.replacement, keep.inbag, sample.fraction,
        no.split.variables, seed, honesty, ci.group.size, split.regularization, alpha)
    
    forest[["ci.group.size"]] <- ci.group.size
    forest[["original.data"]] <- input.data
    forest[["feature.indices"]] <- 1:ncol(X)
    class(forest) <- c("instrumental_forest", "grf")
    forest
}

#' Predict with an instrumental forest
#' 
#' Gets estimates of tau(x) using a trained instrumental forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL,
#'                makes out-of-bag predictions on the training set instead
#'                (i.e., provides predictions at Xi using only trees that did
#'                not use the i-th training example).
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions, along with (optional) variance estimates.
#' @export
predict.instrumental_forest <- function(object, newdata = NULL,
                                        num.threads = NULL, 
                                        estimate.variance = FALSE,
                                        ...) {
    
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
    
    sparse.data <- as.matrix(0)
    variable.names <- character(0)
    
    if (estimate.variance) {
        ci.group.size = object$ci.group.size
    } else {
        ci.group.size = 1
    }
    
    forest.short <- object[-which(names(object) == "original.data")]
    
    if (!is.null(newdata)) {
        input.data <- as.matrix(cbind(newdata, NA))
        instrumental_predict(forest.short, input.data, sparse.data, variable.names, num.threads, 
                             ci.group.size)
    } else {
        input.data <- object[["original.data"]]
        instrumental_predict_oob(forest.short, input.data, sparse.data, variable.names, 
                                 num.threads, ci.group.size)
    }
}
