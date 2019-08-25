validate_X <- function(X) {
  if (inherits(X, "matrix") & !is.numeric(X)) {
    stop(paste(
      "The feature matrix X must be numeric. GRF does not",
      "currently support non-numeric features. If factor variables",
      "are required, we recommend one of the following: Either",
      "represent the factor with a 1-vs-all expansion,",
      "(e.g., using model.matrix(~. , data=X)), or then encode the factor",
      "as a numeric via any natural ordering (e.g., if the factor is a month)."
    ))
  }

  if (inherits(X, "Matrix") & !(inherits(X, "dgCMatrix"))) {
    stop("Currently only sparse data of class 'dgCMatrix' is supported.")
  }

  if (any(is.na(X))) {
    stop("The feature matrix X contains at least one NA.")
  }
}

validate_observations <- function(V, X) {
  if (is.matrix(V) && ncol(V) == 1) {
    V <- as.vector(V)
  } else if (!is.vector(V)) {
    stop(paste("Observations (W, Y, or Z) must be vectors."))
  }

  if (!is.numeric(V) && !is.logical(V)) {
    stop(paste(
      "Observations (W, Y, or Z) must be numeric. GRF does not ",
      "currently support non-numeric observations."
    ))
  }

  if (any(is.na(V))) {
    stop("The vector of observations (W, Y, or Z) contains at least one NA.")
  }

  if (length(V) != nrow(X)) {
    stop("length of observation (W, Y, or Z) does not equal nrow(X).")
  }
  V
}

validate_mtry <- function(mtry, X) {
  if (is.null(mtry)) {
    num.col <- ncol(X)
    default <- ceiling(sqrt(num.col) + 20)
    return(min(default, num.col))
  } else if (!is.numeric(mtry) || mtry < 0) {
    stop("Error: Invalid value for mtry")
  }
  mtry
}

validate_num_threads <- function(num.threads) {
  if (is.null(num.threads)) {
    num.threads <- 0
  } else if (!is.numeric(num.threads) | num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  num.threads
}

validate_min_node_size <- function(min.node.size) {
  if (is.null(min.node.size)) {
    min.node.size <- 5
  } else if (!is.numeric(min.node.size) | min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  min.node.size
}

validate_sample_fraction <- function(sample.fraction) {
  if (is.null(sample.fraction)) {
    sample.fraction <- 0.5
  } else if (!is.numeric(sample.fraction) | sample.fraction <= 0 | sample.fraction > 1) {
    stop("Error: Invalid value for sample.fraction. Please give a value in (0,1].")
  }
  sample.fraction
}

validate_alpha <- function(alpha) {
  if (is.null(alpha)) {
    alpha <- 0.05
  } else if (!is.numeric(alpha) | alpha < 0 | alpha > 0.25) {
    stop("Error: Invalid value for alpha. Please give a value in [0,0.25].")
  }
  alpha
}

validate_imbalance_penalty <- function(imbalance.penalty) {
  if (is.null(imbalance.penalty)) {
    imbalance.penalty <- 0.0
  } else if (!is.numeric(imbalance.penalty) | imbalance.penalty < 0) {
    stop("Error: Invalid value for alpha. Please give a non-negative value.")
  }
  imbalance.penalty
}

validate_seed <- function(seed) {
  if (is.null(seed)) {
    seed <- runif(1, 0, .Machine$integer.max)
  }
  seed
}

validate_clusters <- function(clusters, X) {
  if (is.null(clusters) || length(clusters) == 0) {
    return(vector(mode = "numeric", length = 0))
  }
  if (mode(clusters) != "numeric") {
    stop("Clusters must be able to be coerced to a numeric vector.")
  }
  clusters <- as.numeric(clusters)
  if (!all(clusters == floor(clusters))) {
    stop("Clusters vector cannot contain floating point values.")
  } else if (length(clusters) != nrow(X)) {
    stop("Clusters vector has incorrect length.")
  } else {
    # convert to integers between 0 and n clusters
    clusters <- as.numeric(as.factor(clusters)) - 1
  }
  clusters
}

validate_samples_per_cluster <- function(samples.per.cluster, clusters) {
  if (is.null(clusters) || length(clusters) == 0) {
    return(0)
  }
  cluster_size_counts <- table(clusters)
  min_size <- unname(cluster_size_counts[order(cluster_size_counts)][1])
  if (is.null(samples.per.cluster)) {
    samples.per.cluster <- min_size
  } else if (samples.per.cluster <= 0) {
    stop("samples.per.cluster must be positive")
  }
  samples.per.cluster
}

validate_honesty_fraction <- function(honesty.fraction, honesty) {
  if (!honesty) {
      return(0)
  } else if (is.null(honesty.fraction)) {
    return(0.5)
  } else if (honesty.fraction > 0 && honesty.fraction < 1) {
    return(honesty.fraction)
  } else {
    stop("honesty.fraction must be a positive real number less than 1.")
  }
}

validate_prune_empty_leaves <- function(prune.empty.leaves) {
  if (is.null(prune.empty.leaves)) {
    return(TRUE)
  }
  prune.empty.leaves
}

validate_boost_error_reduction <- function(boost.error.reduction) {
  if (boost.error.reduction < 0 || boost.error.reduction > 1) {
    stop("boost.error.reduction must be between 0 and 1")
  }
  boost.error.reduction
}

validate_ll_vars <- function(linear.correction.variables, num.cols) {
  if (is.null(linear.correction.variables)) {
    linear.correction.variables <- 1:num.cols
  }
  if (min(linear.correction.variables) < 1) {
    stop("Linear correction variables must take positive integer values.")
  } else if (max(linear.correction.variables) > num.cols) {
    stop("Invalid range of correction variables.")
  } else if (!is.vector(linear.correction.variables) |
    !all(linear.correction.variables == floor(linear.correction.variables))) {
    stop("Linear correction variables must be a vector of integers.")
  }
  linear.correction.variables
}

validate_ll_lambda <- function(lambda) {
  if (lambda < 0) {
    stop("Lambda cannot be negative.")
  } else if (!is.numeric(lambda) | length(lambda) > 1) {
    stop("Lambda must be a scalar.")
  }
  lambda
}

validate_ll_path <- function(lambda.path) {
  if (is.null(lambda.path)) {
    lambda.path <- c(0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 1, 10)
  } else if (min(lambda.path) < 0) {
    stop("Lambda values cannot be negative.")
  } else if (!is.numeric(lambda.path)) {
    stop("Lambda values must be numeric.")
  }
  lambda.path
}

validate_newdata <- function(newdata, X) {
  if (ncol(newdata) != ncol(X)) {
    stop("newdata must have the same number of columns as the training matrix.")
  }
  validate_X(newdata)
}

validate_sample_weights <- function(sample.weights, X) {
  if (!is.null(sample.weights)) {
    if (length(sample.weights) != nrow(X)) {
      stop("sample.weights has incorrect length")
    }
    if (any(sample.weights < 0)) {
      stop("sample.weights must be nonnegative")
    }
  }
}

#' @importFrom Matrix Matrix cBind
#' @importFrom methods new
create_data_matrices <- function(X, ..., sample.weights = NULL) {
  default.data <- matrix(nrow = 0, ncol = 0)
  sparse.data <- new("dgCMatrix", Dim = c(0L, 0L))

  if (inherits(X, "dgCMatrix") && ncol(X) > 1) {
    sparse.data <- cbind(X, ..., sample.weights)
  } else {
    X <- as.matrix(X)
    default.data <- as.matrix(cbind(X, ..., sample.weights))
  }

  list(default = default.data, sparse = sparse.data)
}

observation_weights <- function(forest) {
    sample.weights <- if (is.null(forest$sample.weights)) {
        rep(1, length(forest$Y.orig))
    } else {
        forest$sample.weights * length(forest$Y.orig) / sum(forest$sample.weights)
    }
    if (length(forest$clusters) == 0) {
        observation.weight <- sample.weights
    } else {
        clust.factor <- factor(forest$clusters)
        inverse.counts <- 1 / as.numeric(Matrix::colSums(Matrix::sparse.model.matrix(~ clust.factor + 0)))
        observation.weight <- sample.weights * inverse.counts[as.numeric(clust.factor)]
    }
    observation.weight
}
