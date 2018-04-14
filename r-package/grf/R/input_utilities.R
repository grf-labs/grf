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
  # Check for whether this number is too small?
  if (is.null(samples_per_cluster)) {
    samples_per_cluster <- min_size
  } else if (samples_per_cluster > min_size) {
    stop(paste("Smallest cluster has", min_size, "observations",
         "samples_per_cluster of", samples_per_cluster, "is too large."))
  }
  samples_per_cluster
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
