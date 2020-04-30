#' Survival forest
#'
#' Trains a forest for right-censored surival data that can be used to estimate the
#' conditional surival function S(t, x) = P[Y > t | X = x]
#'
#' @param X The covariates.
#' @param Y The response time (may be negative).
#' @param D The censoring indicator (1: failure, 0: censored)
#' @param W The treatment
#' @param num.trees Number of trees grown in the forest. Default is 1000.
#' @param sample.weights (experimental) Weights given to an observation in prediction.
#'                       If NULL, each observation is given the same weight. Default is NULL.
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
#'  parameter tuning, see the grf
#'  \href{https://grf-labs.github.io/grf/REFERENCE.html#honesty-honesty-fraction-honesty-prune-leaves}{algorithm reference}.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. Default is 0.5 (i.e. half of the data is used for
#'                         determining splits).
#' @param honesty.prune.leaves If TRUE, prunes the estimation sample tree such that no leaves
#'  are empty. If FALSE, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to FALSE may improve performance on
#'  small/marginally powered data, but requires more trees (note: tuning does not adjust the number of trees).
#'  Only applies if honesty is enabled. Default is TRUE.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05
#'  (meaning the count of failures on each side of a split has to be at least 5 \% of the total observation count in a node)
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained survival_forest forest object. The attribute `failure.times` contains the unique failure
#'  times in the data set.
#'
#' @examples
#' \dontrun{
#' # Train a standard survival forest.
#' n <- 100
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' failure.time <- -log(runif(n)) * exp(0.1 * X[, 1])
#' censor.time <- rexp(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' s.forest <- survival_forest(X, Y, D)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' s.pred <- predict(s.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' s.pred <- predict(s.forest)
#' }
#'
# ' @export
causal_survival_forest <- function(X, Y, D, W,
                                  tau = max(Y),
                                  num.trees = 1000,
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
  Y <- validate_observations(Y, X)
  D <- validate_observations(D, X)
  if(!all(D %in% c(1, 0))) {
    stop("The censor values can only be 0 or 1.")
  }
  W <- validate_observations(W, X)
  clusters <- validate_clusters(clusters, X)
  samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
  num.threads <- validate_num_threads(num.threads)

  forest.W = regression_forest(X, W)
  W.hat = predict(forest.W)$predictions

  # Estimate m(x) = E[Y|X] = e(X) E[Y | X, W = 1] + (1 - e(x)) E[Y | X, W = 0]
  treated = which(W == 1)
  control = -treated
  sf.treated = survival_forest(X[treated, ], Y[treated], D[treated])
  pred.treated = matrix(0, n, length(sf.treated$failure.times))
  pred.treated[treated, ] = predict(sf.treated)$predictions
  pred.treated[control, ] = predict(sf.treated, X[control, ])$predictions
  expected.treated = expected_survival(pred.treated, sf.treated$failure.times, tau)

  sf.control = survival_forest(X[control, ], Y[control], D[control])
  pred.control = matrix(0, n, length(sf.control$failure.times))
  pred.control[control, ] = predict(sf.control)$predictions
  pred.control[treated, ] = predict(sf.control, X[treated, ])$predictions
  expected.control = expected_survival(pred.control, sf.control$failure.times, tau)

  m = W.hat * expected.treated + (1 - W.hat) * expected.control

  # The rest
  sf.censor = survival_forest(cbind(X, W), Y, 1 - D)
  sf.survival = survival_forest(cbind(X, W), Y, D)

  surv.censor = predict(sf.censor)$predictions
  surv.survival = predict(sf.survival)$predictions
  time.censor = sf.censor$failure.times #"unique_censor"
  time.survival = sf.survival$failure.times
  Y.relabeled.censor = sf.censor$Y.relabeled #"Y_censor_point"

  Y.relabeled.censor[Y.relabeled.censor == 0] = 1 #check
  surv.censor.Y = surv.censor[cbind(1:n, Y.relabeled.censor)]
  surv.censor.Y[sf.censor$Y.relabeled==0] = 1

  START = proc.time()
  # denominator
  denominator = (W - W.hat)^2 * (D == 1) / surv.censor.Y

  Y.relabeled.survival = sf.survival$Y.relabeled
  # extend distribution to all Y points
  Y.relabeled.censor[sf.censor$Y.relabeled==0] = 0 #ch
  surv.censor.extend = cbind(1, surv.censor)[, sort(Y.relabeled.censor) + 1]
  surv.survival.extend = cbind(1, surv.survival)[, sort(Y.relabeled.survival) + 1]

  # numerator
  # Qhat:
  Y.sort = sort(Y)
  N = length(Y.sort)
  at.fail = findInterval(Y, Y.sort)
  q.grid = diff(c(0, Y.sort))

  q.hat = sapply(1:N, function(x) {
    if (at.fail[x] == N || surv.survival.extend[x, at.fail[x]] <=0)
      return (Y.sort[at.fail[x]])
    Y.sort[at.fail[x]] +
      sum(q.grid[(at.fail[x] + 1):N] * surv.survival.extend[x, at.fail[x]:(n - 1)]) /
      surv.survival.extend[x, at.fail[x]] # NB zeros here
  })

  numerator.one = (W - W.hat) *
    ifelse(D == 1, Y - m, q.hat - m) / surv.censor.Y

  # Second numerator
  q.hat.all.t = matrix(NA, N, N) # NxN this will grow very big!!
  swept = sweep(surv.survival.extend[, 1:(N-1), drop = FALSE], 2, q.grid[2:N], "*")
  full.sum = rowSums(swept[, 1:(N-1), drop = FALSE])
  c.expect = full.sum
  q.hat.all.t[, 1] = c.expect
  for (i in 2:(N-1)) {
    c.expect = c.expect - swept[, i-1]
    q.hat.all.t[, i] = c.expect
  }
  q.hat.all.t = q.hat.all.t / surv.survival.extend
  q.hat.all.t[is.na(q.hat.all.t)] = 0
  q.hat.all.t = sweep(q.hat.all.t, 2, Y.sort, "+")
  # as before:
  q.hat.all.t[, N] = max(Y.sort)
  q.hat.all.t = sweep(q.hat.all.t, 1, m, "-") # "q.hat.all.t - as.vector(m)"

  Lambda = t(apply(-log(cbind(1, surv.censor.extend)), 1, diff))
  integrate = sweep(Lambda/surv.censor.extend*q.hat.all.t, 2, q.grid, FUN="*")
  numerator.two = rep(0, N)
  for (i in 1:N) {
    numerator.two[i] = sum(integrate[i, 1:at.fail[i]]) * (W[i] - W.hat[i])
  }

  numerator = numerator.one - numerator.two


  print("Time nuisance: ")
  print(proc.time() - START)
  return (list(numerator=numerator, denominator=denominator, Y=Y, delta=D,surv.censor=surv.censor,
                surv.survival=surv.survival, assign=W, treat=W.hat, expectsurvival1=expected.treated,
              expectsurvival0=expected.control,
            sf.censor=sf.censor,
            sf.survival=sf.survival,
          numerator.one=numerator.one,
        q.hat.all.t=q.hat.all.t,
          numerator.two=numerator.two))
  # return (1)

  forest <- do.call.rcpp(causal_survival_train, c(data, args))
  class(forest) <- c("causal_survival_forest", "grf")
  forest[["X.orig"]] <- X
  forest[["Y.orig"]] <- Y
  forest[["Y.relabeled"]] <- Y.relabeled
  forest[["D.orig"]] <- D #add w.orig
  forest[["sample.weights"]] <- sample.weights
  forest[["clusters"]] <- clusters
  forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
  forest[["has.missing.values"]] <- has.missing.values
  forest[["failure.times"]] <- failure.times

  forest
}

#' Predict with a survival forest forest
#'
#' Gets estimates of the conditional survival function S(t, x) using a trained surival forest (estimated using
#' Kaplan-Meier).
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions.
#'
#' @examples
#' \dontrun{
#' # Train a standard survival forest.
#' n <- 100
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' failure.time <- -log(runif(n)) * exp(0.1 * X[, 1])
#' censor.time <- rexp(n)
#' Y <- pmin(failure.time, censor.time)
#' D <- as.integer(failure.time <= censor.time)
#' s.forest <- survival_forest(X, Y, D)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' s.pred <- predict(s.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' s.pred <- predict(s.forest)
#' }
#'
#' @method predict causal_survival_forest
#' @export
predict.causal_survival_forest <- function(object, newdata = NULL, num.threads = NULL, ...) {
  failure.times = object[["failure.times"]]
  # If possible, use pre-computed predictions.
  if (is.null(newdata) & !is.null(object$predictions)) {
    return(list(predictions = object$predictions, failure.times = failure.times))
  }

  num.threads <- validate_num_threads(num.threads)

  forest.short <- object[-which(names(object) == "X.orig")]
  X <- object[["X.orig"]]
  train.data <- create_train_matrices(X, outcome = object[["Y.relabeled"]], censor = object[["D.orig"]])

  args <- list(forest.object = forest.short,
               num.threads = num.threads,
               num.failures = length(object[["failure.times"]]))

  if (!is.null(newdata)) {
    validate_newdata(newdata, X, allow.na = TRUE)
    test.data <- create_train_matrices(newdata, train.data = FALSE)
    ret <- do.call.rcpp(survival_predict, c(train.data, test.data, args))
  } else {
    ret <- do.call.rcpp(survival_predict_oob, c(train.data, args))
  }

  list(predictions = ret[["predictions"]], failure.times = failure.times)
}

## MISC tmp peripheral

expected_survival = function(predictions, failure.times, tau = NULL) {
  # DOUBLE check these, off by one etc...
  # some facts: E[T] = integral S(t) dt, S(t) = P(T>=t)
  failure.times[failure.times >= tau] = tau

  cbind(1, predictions) %*% diff(c(0, failure.times, tau))
}
