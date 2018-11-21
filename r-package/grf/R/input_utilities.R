validate_X <- function(X) {
  if(inherits(X, "matrix") & !is.numeric(X)) {
    stop(paste("The feature matrix X must numeric. GRF does not", 
         "currently support non-numeric features. If factor variables",
         "are required, we recommend one of the following: Either",
         "represent the factor with a 1-vs-all expansion,",
         "(e.g., using model.matrix(~. , data=X)), or then encode the factor",
         "as a numeric via any natural ordering (e.g., if the factor is a month)."))
  }

  if (inherits(X, "Matrix") & !(inherits(X, "dgCMatrix"))) {
      stop("Currently only sparse data of class 'dgCMatrix' is supported.")
  }
}

validate_mtry <- function(mtry, X) {
  if (is.null(mtry)) {
    num.col = ncol(X)
    default = ceiling(sqrt(num.col) + 20)
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
  if (is.null(clusters)) {
    clusters <- vector(mode="numeric", length=0)
  } else if (length(clusters) == 0) {
    clusters <- vector(mode="numeric", length=0)
  } else if (!is.vector(clusters) | !all(clusters == floor(clusters))) {
    stop("Clusters must be a vector of integers.")
  } else if (length(clusters) != nrow(X)) {
    stop("Clusters has incorrect length.")
  } else {
    # convert to integers between 0 and n clusters
    clusters <- as.numeric(as.factor(clusters)) - 1
  }
  clusters
}

validate_samples_per_cluster <- function(samples_per_cluster, clusters) {
  if (is.null(clusters) || length(clusters) == 0) {
    return(0)
  }
  cluster_size_counts <- table(clusters)
  min_size <- unname(cluster_size_counts[order(cluster_size_counts)][1])
  if (is.null(samples_per_cluster)) {
    samples_per_cluster <- min_size
  } else if (samples_per_cluster <= 0) {
    stop("samples_per_cluster must be positive")
  }
  samples_per_cluster
}

validate_honesty_fraction <- function(honesty.fraction, honesty) {
  if (!honesty) {
      if (is.null(honesty.fraction)) {
        return(NULL)
      }
      else {
        stop("honesty.fraction is not used when honesty = FALSE and should be NULL in this case.")
      }
  } else if (is.null(honesty.fraction)) {
    return(0.5)
  } else if (honesty.fraction > 0 && honesty.fraction < 1) {
    return(honesty.fraction)
  } else {
    stop("honesty.fraction must be a positive real number less than 1.")
  }
}

validate_ll_vars <- function(linear.correction.variables, num.cols){
  if (is.null(linear.correction.variables)) {
    linear.correction.variables = 1:num.cols
  }
  if (min(linear.correction.variables) < 1) {
    stop("Linear correction variables must take positive integer values.")
  } else if (max(linear.correction.variables) > num.cols) {
    stop("Invalid range of correction variables.")
  } else if (!is.vector(linear.correction.variables) | !all(linear.correction.variables == floor(linear.correction.variables))) {
    stop("Linear correction variables must be a vector of integers.")
  }
  linear.correction.variables
}

validate_ll_lambda <- function(lambda){
  if (lambda < 0) {
    stop("Lambda cannot be negative.")
  } else if (!is.numeric(lambda) | length(lambda) > 1) {
    stop("Lambda must be a scalar.")
  }
  lambda
}

validate_ll_path <- function(lambda.path){
  if (is.null(lambda.path)) {
    lambda.path = c(0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 1, 10)
  } else if (min(lambda.path)<0) {
    stop("Lambda values cannot be negative.")
  } else if (!is.numeric(lambda.path)) {
    stop("Lambda values must be numeric.")
  }
  lambda.path
}

coerce_honesty_fraction <- function(honesty.fraction) {
  if(is.null(honesty.fraction)) {
    return(0)
  }
  honesty.fraction
}

create_data_matrices <- function(X, ...) {
  default.data <- matrix(nrow=0, ncol=0);    
  sparse.data <- new("dgCMatrix", Dim = c(0L, 0L))

  if (inherits(X, "dgCMatrix") && ncol(X) > 1) {
    sparse.data <- cbind(X, ...)
  } else {
    default.data <- as.matrix(cbind(X, ...))
  }

  list(default = default.data, sparse = sparse.data)
}
