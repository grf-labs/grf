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
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. If Y.hat = NULL, these are estimated using
#'              a separate regression forest.
#' @param W.hat Estimates of the treatment propensities E[W | Xi]. If W.hat = NULL,
#'              these are estimated using a separate regression forest.
#' @param Z.hat Estimates of the instrument propensities E[Z | Xi]. If Z.hat = NULL,
#'              these are estimated using a separate regression forest.
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
#' @param ci.group.size The forst will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2.
#' @param reduced.form.weight Whether splits should be regularized towards a naive
#'                            splitting criterion that ignores the instrument (and
#'                            instead emulates a causal forest).
#' @param alpha A tuning parameter that controls the maximum imbalance of a split.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized.
#' @param stabilize.splits Whether or not the instrument should be taken into account when
#'                         determining the imbalance of a split (experimental).
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed.
#' @param seed The seed for the C++ random number generator.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#' @param samples_per_cluster If sampling by cluster, the number of observations to be sampled from
#'                            each cluster when training a tree. If NULL, we set samples_per_cluster to the size
#'                            of the smallest cluster. If some clusters are smaller than samples_per_cluster,
#'                            the whole cluster is used every time the cluster is drawn. Note that
#'                            clusters with less than samples_per_cluster observations get relatively
#'                            smaller weight than others in training the forest, i.e., the contribution
#'                            of a given cluster to the final forest scales with the minimum of
#'                            the number of observations in the cluster and samples_per_cluster.
#'
#' @return A trained instrumental forest object.
#' @export
instrumental_forest <- function(X, Y, W, Z,
                                Y.hat = NULL,
                                W.hat = NULL,
                                Z.hat = NULL,
                                sample.fraction = 0.5,
                                mtry = NULL,
                                num.trees = 2000,
                                num.threads = NULL,
                                min.node.size = NULL,
                                honesty = TRUE,
                                honesty.fraction = NULL,
                                ci.group.size = 2,
                                reduced.form.weight = 0,
                                alpha = 0.05,
                                imbalance.penalty = 0.0,
                                stabilize.splits = TRUE,
                                compute.oob.predictions = TRUE,
                                seed = NULL,
                                clusters = NULL,
                                samples_per_cluster = NULL) {
    validate_X(X)
    if(length(Y) != nrow(X)) { stop("Y has incorrect length.") }
    if(length(W) != nrow(X)) { stop("W has incorrect length.") }
    if(length(Z) != nrow(X)) { stop("Z has incorrect length.") }
    
    
    mtry <- validate_mtry(mtry, X)
    num.threads <- validate_num_threads(num.threads)
    min.node.size <- validate_min_node_size(min.node.size)
    sample.fraction <- validate_sample_fraction(sample.fraction)
    seed <- validate_seed(seed)
    clusters <- validate_clusters(clusters, X)
    samples_per_cluster <- validate_samples_per_cluster(samples_per_cluster, clusters)
    honesty.fraction <- validate_honesty_fraction(honesty.fraction, honesty)
    
    if (!is.numeric(reduced.form.weight) | reduced.form.weight < 0 | reduced.form.weight > 1) {
        stop("Error: Invalid value for reduced.form.weight. Please give a value in [0,1].")
    }
    
    if (is.null(Y.hat)) {
      forest.Y <- regression_forest(X, Y, sample.fraction = sample.fraction, mtry = mtry, 
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      Y.hat <- predict(forest.Y)$predictions
    } else if (length(Y.hat) == 1) {
      Y.hat <- rep(Y.hat, nrow(X))
    } else if (length(Y.hat) != nrow(X)) {
      stop("Y.hat has incorrect length.")
    }
    
    if (is.null(W.hat)) {
      forest.W <- regression_forest(X, W, sample.fraction = sample.fraction, mtry = mtry, 
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      W.hat <- predict(forest.W)$predictions
    } else if (length(W.hat) == 1) {
      W.hat <- rep(W.hat, nrow(X))
    } else if (length(W.hat) != nrow(X)) {
      stop("W.hat has incorrect length.")
    }
    
    if (is.null(Z.hat)) {
      forest.Z <- regression_forest(X, Z, sample.fraction = sample.fraction, mtry = mtry, 
                                    num.trees = min(500, num.trees), num.threads = num.threads, min.node.size = NULL, honesty = TRUE,
                                    honesty.fraction = NULL, seed = seed, ci.group.size = 1, alpha = alpha, imbalance.penalty = imbalance.penalty,
                                    clusters = clusters, samples_per_cluster = samples_per_cluster);
      Z.hat <- predict(forest.Z)$predictions
    } else if (length(Z.hat) == 1) {
      Z.hat <- rep(Z.hat, nrow(X))
    } else if (length(Z.hat) != nrow(X)) {
      stop("Z.hat has incorrect length.")
    }
    
    data <- create_data_matrices(X, Y - Y.hat, W - W.hat, Z - Z.hat)
    
    outcome.index <- ncol(X) + 1
    treatment.index <- ncol(X) + 2
    instrument.index <- ncol(X) + 3
    
    forest <- instrumental_train(data$default, data$sparse, outcome.index, treatment.index,
        instrument.index, mtry, num.trees, num.threads, min.node.size, sample.fraction, seed, honesty,
        coerce_honesty_fraction(honesty.fraction), ci.group.size, reduced.form.weight, alpha, 
        imbalance.penalty, stabilize.splits, clusters, samples_per_cluster)

    forest[["ci.group.size"]] <- ci.group.size
    forest[["X.orig"]] <- X
    forest[["Y.orig"]] <- Y
    forest[["W.orig"]] <- W
    forest[["Z.orig"]] <- Z
    forest[["Y.hat"]] <- Y.hat
    forest[["W.hat"]] <- W.hat
    forest[["Z.hat"]] <- Z.hat
    forest[["clusters"]] <- clusters

    class(forest) <- c("instrumental_forest", "grf")

    if (compute.oob.predictions) {
        oob.pred <- predict(forest)
        forest[["predictions"]] <- oob.pred$predictions
        forest[["debiased.error"]] <- oob.pred$debiased.error
    }

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
#'
#' @method predict instrumental_forest
#' @export
predict.instrumental_forest <- function(object, newdata = NULL,
                                        num.threads = NULL, 
                                        estimate.variance = FALSE,
                                        ...) {

    # If possible, use pre-computed predictions.
    if (is.null(newdata) & !estimate.variance & !is.null(object$predictions)) {
        return(data.frame(predictions=object$predictions,
                          debiased.error=object$debiased.error))
    }

    num.threads <- validate_num_threads(num.threads)

    if (estimate.variance) {
        ci.group.size <- object$ci.group.size
    } else {
        ci.group.size <- 1
    }

    forest.short <- object[-which(names(object) == "X.orig")]

    if (!is.null(newdata)) {
        data <- create_data_matrices(newdata)
        ret <- instrumental_predict(forest.short, data$default, data$sparse,
                                    num.threads, ci.group.size)
    } else {
        data <- create_data_matrices(object[["X.orig"]])
        ret <- instrumental_predict_oob(forest.short, data$default, data$sparse,
                                        num.threads, ci.group.size)
    }

    # Convert list to data frame.
    empty = sapply(ret, function(elem) length(elem) == 0)
    do.call(cbind.data.frame, ret[!empty])
}
