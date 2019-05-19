get_initial_params <- function(min.node.size,
                               sample.fraction,
                               mtry,
                               alpha,
                               imbalance.penalty) {
 c(min.node.size = if (is.null(min.node.size)) NA else validate_min_node_size(min.node.size),
    sample.fraction = if (is.null(sample.fraction)) NA else validate_sample_fraction(sample.fraction),
    mtry = if (is.null(mtry)) NA else validate_mtry(mtry),
    alpha = if (is.null(alpha)) NA else validate_alpha(alpha),
    imbalance.penalty = if (is.null(imbalance.penalty)) NA else validate_imbalance_penalty(imbalance.penalty))
}

get_params_from_draw <- function(X, draws) {
  if (is.vector(draws))
    draws = rbind(c(draws))
  n = nrow(draws)
  vapply(colnames(draws), function(param) {
    if (param == "min.node.size")
      return(floor(2^(draws[, param] * (log(nrow(X)) / log(2) - 4))))
    else if (param == "sample.fraction")
      return(0.05 + 0.45 * draws[, param])
    else if (param == "mtry")
      return(ceiling(min(ncol(X), sqrt(ncol(X)) + 20) * draws[, param]))
    else if (param == "alpha")
      return(draws[, param] / 4)
    else if (param == "imbalance.penalty")
      return (-log(draws[, param]))
    else
      stop("Unrecognized parameter name provided: ", param)
  }, FUN.VALUE = numeric(n))
}
