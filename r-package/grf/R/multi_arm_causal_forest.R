#' Multi-arm causal forest
#'
#' Trains a causal forest that can be used to estimate
#' conditional average treatment effects tau_k(X). When
#' the treatment assignment W is \{1, ..., K\} and unconfounded,
#' we have tau_k(X) = E[Y(k) - Y(1) | X = x] where Y(k) and
#' Y(1) are potential outcomes corresponding to the treatment
#' state for arm k and the baseline arm 1.
#'
#'
#' This forest fits a multi-arm treatment estimate following the multivariate
#' extension of the "R-learner" suggested in Nie and Wager (2021), with kernel
#' weights derived by the GRF algortim (Athey, Tibshirani, and Wager, 2019).
#' In particular, with K arms, and W encoded as \{0, 1\}^(K-1), we estimate, for
#' a target sample x, and a chosen baseline arm:
#'
#' \eqn{\hat \tau(x) = argmin_{\tau} \left\{ \sum_{i=1}^{n}
#' \alpha_i (x) \left( Y_i - \hat m^{(-i)}(X_i) - c(x) -
#' \left\langle W_i - \hat e^{(-i)}(X_i), \,  \tau(X_i)  \right\rangle
#' \right)^2 \right\}},
#'
#' where the angle brackets indicates an inner product, e(X) = E[W | X = x] is
#' a (vector valued) generalized propensity score, and m(x) = E[Y | X = x].
#' The forest weights alpha(x) are derived from a generalized random forest
#' splitting on the vector-valued gradient of tau(x). (The intercept c(x)
#' is a nuisance parameter not directly estimated). By default, e(X) and
#' m(X) are estimated using two separate random forests, a probability forest
#' and regression forest respectively (optionally provided through the arguments
#' W.hat and Y.hat). The k-th element of tau(x) measures the conditional average
#' treatment effect of the k-th treatment arm at X = x for k = 1, ..., K-1.
#' The treatment effect for multiple outcomes can be estimated jointly (i.e. Y can
#' be vector-valued) - in which case the splitting rule takes into account
#' all outcomes simultaneously (specifically, we concatenate the gradient
#' vector for each outcome).
#'
#' For a single treatment, this forest is equivalent to a causal forest, however,
#' they may produce different results due to differences in numerics.
#'
#' @param X The covariates used in the causal regression.
#' @param Y The outcome (must be a numeric vector or matrix [one column per outcome] with no NAs).
#'  Multiple outcomes should be on the same scale.
#' @param W The treatment assignment (must be a factor vector with no NAs). The reference treatment
#'          is set to the first treatment according to the ordinality of the factors, this can be changed
#'          with the `relevel` function in R.
#' @param Y.hat Estimates of the expected responses E[Y | Xi], marginalizing
#'              over treatment. If Y.hat = NULL, these are estimated using
#'              a separate multi-task regression forest. Default is NULL.
#' @param W.hat Matrix with estimates of the treatment propensities E[Wk | Xi].
#'              If W.hat = NULL, these are estimated using a probability forest.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to each sample in estimation.
#'                       If NULL, each observation receives the same weight.
#'                       Note: To avoid introducing confounding, weights should be
#'                       independent of the potential outcomes given X. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param equalize.cluster.weights If FALSE, each unit is given the same weight (so that bigger
#'  clusters get more weight). If TRUE, each cluster is given equal weight in the forest. In this case,
#'  during training, each tree uses the same number of observations from each drawn cluster: If the
#'  smallest cluster has K units, then when we sample a cluster during training, we only give a random
#'  K elements of the cluster to the tree-growing procedure. When estimating average treatment effects,
#'  each observation is given weight 1/cluster size, so that the total weight of each cluster is the
#'  same. Note that, if this argument is FALSE, sample weights may also be directly adjusted via the
#'  sample.weights argument. If this argument is TRUE, sample.weights must be set to NULL. Default is
#'  FALSE.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction. Default is 0.5.
#' @param mtry Number of variables tried for each split. Default is
#'             \eqn{\sqrt p + 20} where p is the number of variables.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#'                      Default is 5.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting). Default is TRUE.
#'  For a detailed description of honesty, honesty.fraction, honesty.prune.leaves, and recommendations for
#'  parameter tuning, see the grf algorithm reference.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. Default is 0.5 (i.e. half of the data is used for
#'                         determining splits).
#' @param honesty.prune.leaves If TRUE, prunes the estimation sample tree such that no leaves
#'  are empty. If FALSE, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to FALSE may improve performance on
#'  small/marginally powered data, but requires more trees (note: tuning does not adjust the number of trees).
#'  Only applies if honesty is enabled. Default is TRUE.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized. Default is 0.
#' @param stabilize.splits Whether or not the treatment should be taken into account when
#' determining the imbalance of a split. It is an exact extension of the single-arm constraints (detailed
#' in the causal forest algorithm reference) to multiple arms, where the constraints apply to each treatment arm independently.
#' Default is TRUE.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2. (Confidence intervals are
#'                      currently only supported for univariate outcomes Y).
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained multi arm causal forest object.
#'
#' @references Athey, Susan, Julie Tibshirani, and Stefan Wager. "Generalized Random Forests".
#'  Annals of Statistics, 47(2), 2019.
#' @references Nie, Xinkun, and Stefan Wager. "Quasi-Oracle Estimation of Heterogeneous Treatment Effects".
#'  Biometrika, 108(2), 2021.
#'
#' @examples
#' \donttest{
#' # Train a multi arm causal forest.
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' Y <- X[, 1] + X[, 2] * (W == "B") - 1.5 * X[, 2] * (W == "C") + rnorm(n)
#' mc.forest <- multi_arm_causal_forest(X, Y, W)
#'
#' # Predict contrasts (out-of-bag) using the forest.
#' # By default, the first ordinal treatment is used as baseline ("A" in this example),
#' # giving two contrasts tau_B = Y(B) - Y(A), tau_C = Y(C) - Y(A)
#' mc.pred <- predict(mc.forest)
#' # Fitting several outcomes jointly is supported, and the returned prediction array has
#' # dimension [num.samples, num.contrasts, num.outcomes]. Since num.outcomes is one in
#' # this example, we can `drop()` this singleton dimension using the `[,,]` shorthand.
#' tau.hat <- mc.pred$predictions[,,]
#'
#' plot(X[, 2], tau.hat[, "B - A"], ylab = "tau.contrast")
#' abline(0, 1, col = "red")
#' points(X[, 2], tau.hat[, "C - A"], col = "blue")
#' abline(0, -1.5, col = "red")
#' legend("topleft", c("B - A", "C - A"), col = c("black", "blue"), pch = 19)
#'
#' # The average treatment effect of the arms with "A" as baseline.
#' average_treatment_effect(mc.forest)
#'
#' # The conditional response surfaces mu_k(X) for a single outcome can be reconstructed from
#' # the contrasts tau_k(x), the treatment propensities e_k(x), and the conditional mean m(x).
#' # Given treatment "A" as baseline we have:
#' # m(x) := E[Y | X] = E[Y(A) | X] + E[W_B (Y(B) - Y(A))] + E[W_C (Y(C) - Y(A))]
#' # which given unconfoundedness is equal to:
#' # m(x) = mu(A, x) + e_B(x) tau_B(X) + e_C(x) tau_C(x)
#' # Rearranging and plugging in the above expressions, we obtain the following estimates
#' # * mu(A, x) = m(x) - e_B(x) tau_B(x) - e_C(x) tau_C(x)
#' # * mu(B, x) = m(x) + (1 - e_B(x)) tau_B(x) - e_C(x) tau_C(x)
#' # * mu(C, x) = m(x) - e_B(x) tau_B(x) + (1 - e_C(x)) tau_C(x)
#' Y.hat <- mc.forest$Y.hat
#' W.hat <- mc.forest$W.hat
#'
#' muA <- Y.hat - W.hat[, "B"] * tau.hat[, "B - A"] - W.hat[, "C"] * tau.hat[, "C - A"]
#' muB <- Y.hat + (1 - W.hat[, "B"]) * tau.hat[, "B - A"] - W.hat[, "C"] * tau.hat[, "C - A"]
#' muC <- Y.hat - W.hat[, "B"] * tau.hat[, "B - A"] + (1 - W.hat[, "C"]) * tau.hat[, "C - A"]
#'
#' # These can also be obtained with some array manipulations.
#' # (the first column is always the baseline arm)
#' Y.hat.baseline <- Y.hat - rowSums(W.hat[, -1, drop = FALSE] * tau.hat)
#' mu.hat.matrix <- cbind(Y.hat.baseline, c(Y.hat.baseline) + tau.hat)
#' colnames(mu.hat.matrix) <- levels(W)
#' head(mu.hat.matrix)
#'
#' # The reference level for contrast prediction can be changed with `relevel`.
#' # Fit and predict with treatment B as baseline:
#' W <- relevel(W, ref = "B")
#' mc.forest.B <- multi_arm_causal_forest(X, Y, W)
#' }
#'
#' @export
multi_arm_causal_forest <- function(X, Y, W,
                                    Y.hat = NULL,
                                    W.hat = NULL,
                                    num.trees = 2000,
                                    sample.weights = NULL,
                                    clusters = NULL,
                                    equalize.cluster.weights = FALSE,
                                    sample.fraction = 0.5,
                                    mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                                    min.node.size = 5,
                                    honesty = TRUE,
                                    honesty.fraction = 0.5,
                                    honesty.prune.leaves = TRUE,
                                    alpha = 0.05,
                                    imbalance.penalty = 0,
                                    stabilize.splits = TRUE,
                                    ci.group.size = 2,
                                    compute.oob.predictions = TRUE,
                                    num.threads = NULL,
                                    seed = runif(1, 0, .Machine$integer.max)) {
  has.missing.values <- validate_X(X, allow.na = TRUE)
  validate_sample_weights(sample.weights, X)
  Y <- validate_observations(Y, X, allow.matrix = TRUE)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)
  if (length(W) != nrow(X)) {
    stop("length of observation (W, Y, Z or D) does not equal nrow(X).")
  }
  if (anyNA(W)) {
    stop("The vector of observations (W, Y, Z or D) contains at least one NA.")
  }
  if (!is.factor(W)) {
    stop("The treatment assignment W must be a factor vector.")
  }
  if (nlevels(W) == 1) {
    stop("Can not compute contrasts from a single treatment.")
  }
  if (nlevels(W) != nlevels(droplevels(W))) {
    warning("The treatment vector W contains unused levels (see `droplevels()` to drop unused levels).")
  }

  args.orthog <- list(X = X,
                      num.trees = max(50, num.trees / 4),
                      sample.weights = sample.weights,
                      clusters = clusters,
                      equalize.cluster.weights = equalize.cluster.weights,
                      sample.fraction = sample.fraction,
                      mtry = mtry,
                      min.node.size = 5,
                      honesty = TRUE,
                      honesty.fraction = 0.5,
                      honesty.prune.leaves = honesty.prune.leaves,
                      alpha = alpha,
                      imbalance.penalty = imbalance.penalty,
                      num.threads = num.threads,
                      seed = seed)

  if (is.null(Y.hat)) {
    forest.Y <- do.call(multi_regression_forest, c(Y = list(Y), args.orthog))
    Y.hat <- predict(forest.Y)$predictions
  } else if (is.numeric(Y.hat) && length(Y.hat) == NCOL(Y)) {
    Y.hat <- matrix(Y.hat, nrow = NROW(Y), ncol = NCOL(Y), byrow = TRUE)
  } else if (!(is.matrix(Y.hat) || is.data.frame(Y.hat) || is.vector(Y.hat))) {
      stop("Y.hat should be a matrix of E[Y | Xi] estimates.")
  } else {
    Y.hat <- as.matrix(Y.hat)
    if (NROW(Y.hat) != nrow(X) || NCOL(Y.hat) != NCOL(Y)) {
      stop("Y.hat has incorrect dimensions.")
    }
  }

  if (is.null(W.hat)) {
    args.orthog$ci.group.size <- 1
    forest.W <- do.call(probability_forest, c(Y = list(W), args.orthog))
    W.hat <- predict(forest.W)$predictions
  } else if (is.numeric(W.hat) && length(W.hat) == nlevels(W)) {
    W.hat <- matrix(W.hat, nrow = length(W), ncol = nlevels(W), byrow = TRUE)
  } else if (!(is.matrix(W.hat) || is.data.frame(W.hat))) {
    stop("W.hat should a matrix of E[W_k | Xi] estimates.")
  } else {
    W.hat <- as.matrix(W.hat)
    if ((NROW(W.hat) != nrow(X)) || NCOL(W.hat) != nlevels(W)) {
      stop("W.hat has incorrect dimensions: should be a matrix of E[W_k | Xi] estimates.")
    }
    if (!is.null(colnames(W.hat))) {
      if (!identical(levels(W), colnames(W.hat))) {
        warning(paste(
          "Column names are provided for W.hat, but do not correspond to the treatment levels W",
          "(multi_arm_causal_forest assume the following is TRUE: `identical(levels(W), colnames(W.hat))`)."))
      }
    }
  }

  W.matrix <- stats::model.matrix(~ W - 1)
  Y.centered <- Y - Y.hat
  W.centered <- W.matrix - W.hat
  data <- create_train_matrices(X,
                                outcome = Y.centered,
                                # We always use the first treatment as the reference level
                                treatment = W.centered[, -1],
                                sample.weights = sample.weights)
  args <- list(num.trees = num.trees,
               clusters = clusters,
               samples.per.cluster = samples.per.cluster,
               sample.fraction = sample.fraction,
               mtry = mtry,
               min.node.size = min.node.size,
               honesty = honesty,
               honesty.fraction = honesty.fraction,
               honesty.prune.leaves = honesty.prune.leaves,
               alpha = alpha,
               imbalance.penalty = imbalance.penalty,
               stabilize.splits = stabilize.splits,
               ci.group.size = ci.group.size,
               compute.oob.predictions = compute.oob.predictions,
               num.threads = num.threads,
               seed = seed)

  forest <- do.call.rcpp(multi_causal_train, c(data, args))
  class(forest) <- c("multi_arm_causal_forest", "grf")
  forest[["seed"]] <- seed
  forest[["ci.group.size"]] <- ci.group.size
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["W.orig"]] <- W
  forest[["Y.hat"]] <- Y.hat
  forest[["W.hat"]] <- W.hat
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["sample.weights"]] <- sample.weights
  forest[["has.missing.values"]] <- has.missing.values

  forest
}

#' Predict with a multi arm causal forest
#'
#' Gets estimates of contrasts tau_k(x) using a trained multi arm causal forest (k = 1,...,K-1
#' where K is the number of treatments).
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals). This option is currently
#'                          only supported for univariate outcomes Y.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list with elements `predictions`: a 3d array of dimension [num.samples, K-1, M] with
#' predictions for each contrast, for each outcome 1,..,M (singleton dimensions in this array can
#' be dropped by passing the `drop` argument to `[`, or with the shorthand `$predictions[,,]`),
#'  and optionally `variance.estimates`: a matrix with K-1 columns with variance estimates for each contrast.
#'
#' @examples
#' \donttest{
#' # Train a multi arm causal forest.
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' Y <- X[, 1] + X[, 2] * (W == "B") - 1.5 * X[, 2] * (W == "C") + rnorm(n)
#' mc.forest <- multi_arm_causal_forest(X, Y, W)
#'
#' # Predict contrasts (out-of-bag) using the forest.
#' # By default, the first ordinal treatment is used as baseline ("A" in this example),
#' # giving two contrasts tau_B = Y(B) - Y(A), tau_C = Y(C) - Y(A)
#' mc.pred <- predict(mc.forest)
#' # Fitting several outcomes jointly is supported, and the returned prediction array has
#' # dimension [num.samples, num.contrasts, num.outcomes]. Since num.outcomes is one in
#' # this example, we can `drop()` this singleton dimension using the `[,,]` shorthand.
#' tau.hat <- mc.pred$predictions[,,]
#'
#' plot(X[, 2], tau.hat[, "B - A"], ylab = "tau.contrast")
#' abline(0, 1, col = "red")
#' points(X[, 2], tau.hat[, "C - A"], col = "blue")
#' abline(0, -1.5, col = "red")
#' legend("topleft", c("B - A", "C - A"), col = c("black", "blue"), pch = 19)
#'
#' # The average treatment effect of the arms with "A" as baseline.
#' average_treatment_effect(mc.forest)
#'
#' # The conditional response surfaces mu_k(X) for a single outcome can be reconstructed from
#' # the contrasts tau_k(x), the treatment propensities e_k(x), and the conditional mean m(x).
#' # Given treatment "A" as baseline we have:
#' # m(x) := E[Y | X] = E[Y(A) | X] + E[W_B (Y(B) - Y(A))] + E[W_C (Y(C) - Y(A))]
#' # which given unconfoundedness is equal to:
#' # m(x) = mu(A, x) + e_B(x) tau_B(X) + e_C(x) tau_C(x)
#' # Rearranging and plugging in the above expressions, we obtain the following estimates
#' # * mu(A, x) = m(x) - e_B(x) tau_B(x) - e_C(x) tau_C(x)
#' # * mu(B, x) = m(x) + (1 - e_B(x)) tau_B(x) - e_C(x) tau_C(x)
#' # * mu(C, x) = m(x) - e_B(x) tau_B(x) + (1 - e_C(x)) tau_C(x)
#' Y.hat <- mc.forest$Y.hat
#' W.hat <- mc.forest$W.hat
#'
#' muA <- Y.hat - W.hat[, "B"] * tau.hat[, "B - A"] - W.hat[, "C"] * tau.hat[, "C - A"]
#' muB <- Y.hat + (1 - W.hat[, "B"]) * tau.hat[, "B - A"] - W.hat[, "C"] * tau.hat[, "C - A"]
#' muC <- Y.hat - W.hat[, "B"] * tau.hat[, "B - A"] + (1 - W.hat[, "C"]) * tau.hat[, "C - A"]
#'
#' # These can also be obtained with some array manipulations.
#' # (the first column is always the baseline arm)
#' Y.hat.baseline <- Y.hat - rowSums(W.hat[, -1, drop = FALSE] * tau.hat)
#' mu.hat.matrix <- cbind(Y.hat.baseline, c(Y.hat.baseline) + tau.hat)
#' colnames(mu.hat.matrix) <- levels(W)
#' head(mu.hat.matrix)
#'
#' # The reference level for contrast prediction can be changed with `relevel`.
#' # Fit and predict with treatment B as baseline:
#' W <- relevel(W, ref = "B")
#' mc.forest.B <- multi_arm_causal_forest(X, Y, W)
#' }
#'
#' @method predict multi_arm_causal_forest
#' @export
predict.multi_arm_causal_forest <- function(object,
                                            newdata = NULL,
                                            num.threads = NULL,
                                            estimate.variance = FALSE,
                                            ...) {
  if (estimate.variance && NCOL(object[["Y.orig"]]) > 1) {
    stop("Pointwise variance estimates are only supported for one outcome.")
  }
  treatment.names <- levels(object[["W.orig"]])
  contrast.names <- paste(treatment.names[-1], "-", treatment.names[1])
  outcome.names <- if (is.null(colnames(object[["Y.orig"]]))) {
    paste("Y", 1:NCOL(object[["Y.orig"]]), sep = ".")
  } else {
    make.names(colnames(object[["Y.orig"]]), unique = TRUE)
  }
  # Note the term `num.treatments` is overloaded, in multi_arm_causal_forest's context it means `num.contrasts`
  num.treatments <- length(treatment.names) - 1
  num.outcomes <- length(outcome.names)
  dimnames <- list(NULL, contrast.names, outcome.names)
  # If possible, use pre-computed predictions.
  if (is.null(newdata) && !estimate.variance && !is.null(object$predictions)) {
    predictions <- array(object$predictions, dim = c(NROW(object$predictions), num.treatments, num.outcomes),
                         dimnames = dimnames)
    return(list(predictions = predictions))
  }

  num.threads <- validate_num_threads(num.threads)
  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X)

  args <- list(forest.object = forest.short,
               num.outcomes = num.outcomes,
               num.treatments = num.treatments,
               num.threads = num.threads,
               estimate.variance = estimate.variance)

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_test_matrices(newdata)
    ret <- do.call.rcpp(multi_causal_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(multi_causal_predict_oob, c(train.data, args))
  }
  predictions <- array(ret$predictions, dim = c(NROW(ret$predictions), num.treatments, num.outcomes),
                       dimnames = dimnames)

  list(predictions = predictions,
       variance.estimates = if (estimate.variance) ret$variance.estimates)
}
