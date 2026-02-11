#' grf package options
#'
#' grf package options can be set using R's \code{\link{options}} command.
#' The current available options are:
#' \itemize{
#' \item `grf.progress.bar`: controls whether a progress bar is displayed when training a forest.
#'  The default value is `FALSE`.
#'  \item `grf.legacy.seed`: controls whether grf's random seed behavior depends on
#'  the number of CPU threads used to train the forest. The default value is `FALSE`.
#'  Set to `TRUE` to recover results produced with grf versions prior to 2.4.0.
#' }
#'
#' @return Prints the current grf package options.
#'
#' @examples
#' \donttest{
#' # Enable progress bar when training forests.
#' options(grf.progress.bar = TRUE)
#' n <- 1500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' r.forest <- regression_forest(X, Y)
#'
#' # Print current package options.
#' options(grf.progress.bar = FALSE)
#' grf_options()
#'
#' # Use random seed behavior prior to version 2.4.0.
#' options(grf.legacy.seed = TRUE)
#'
#' # Use random seed independent of num.threads (default as of version 2.4.0 and higher).
#' options(grf.legacy.seed = FALSE)
#' }
#'
#' @export
grf_options <- function() {
    print(c(
        grf.progress.bar = get_progress_bar(),
        grf.legacy.seed = get_legacy_seed()
    ))
}
