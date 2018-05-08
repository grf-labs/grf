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

get_params_from_draw <- function(X, draw) {
  result = c()
  for (param in names(draw)) {
    if (param == "min.node.size") {
      value = floor(2^(draw[param] * (log(nrow(X)) / log(2) - 4)))
    } else if (param == "sample.fraction") {
      value = 0.05 + 0.45 * draw[param]
    } else if (param == "mtry") {
      value = ceiling(min(ncol(X), sqrt(ncol(X)) + 20) * draw[param])
    } else if (param == "alpha") {
      value = draw[param]/4
    } else if (param == "imbalance.penalty") {
      value = -log(draw[param])
    } else {
      stop("Unrecognized parameter name provided: ", param)
    }
    result = c(result, value)
  }
  result
}
