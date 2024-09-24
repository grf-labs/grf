#' grf package options
#'
#' grf package options can be set using R's \code{\link{options}} command.
#' The current available options are:
#' \itemize{
#'  \item `grf.legacy.seed`: controls whether grf's random seed behavior depends on
#'  the number of CPU threads used to train the forest. The default value is `FALSE`.
#'  Set to `TRUE` to recover results produced with grf versions prior to 2.4.0.
#' }
#'
#' @return Prints the current grf package options.
#'
#' @examples
#' \donttest{
#' # Use random seed behavior prior to version 2.4.0.
#' options(grf.legacy.seed = TRUE)
#'
#' # Print current package options.
#' grf_options()
#'
#' # Use random seed independent of num.threads (default as of version 2.4.0 and higher).
#' options(grf.legacy.seed = FALSE)
#' }
#'
#' @export
grf_options <- function() {
    print(c(
        grf.legacy.seed = get_legacy_seed()
    ))
}
