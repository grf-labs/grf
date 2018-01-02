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
    min.node.size <- 0
  } else if (!is.numeric(min.node.size) | min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  min.node.size
}

validate_sample_fraction <- function(sample.fraction) {
  if (!is.numeric(sample.fraction) | sample.fraction <= 0 | sample.fraction > 1) {
    stop("Error: Invalid value for sample.fraction. Please give a value in (0,1].")
  }
  sample.fraction
}

validate_seed <- function(seed) {
  if (is.null(seed)) {
    seed <- runif(1, 0, .Machine$integer.max)
  }
  seed
}

create_data_matrices <- function(X, ...) {
  default.data <- matrix(nrow=0, ncol=0);    
  sparse.data <- new("dgCMatrix", Dim = c(0L, 0L))
  if (inherits(X, "dgCMatrix")) {
    sparse.data = Matrix::cBind(X, ...)
  } else {
    default.data <- as.matrix(cbind(X, ...))
  }

  list(default = default.data, sparse = sparse.data)
}
