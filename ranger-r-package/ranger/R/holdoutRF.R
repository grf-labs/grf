# -------------------------------------------------------------------------------
#   This file is part of Ranger.
#
# Ranger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ranger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ranger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Marvin N. Wright
# Institut fuer Medizinische Biometrie und Statistik
# Universitaet zu Luebeck
# Ratzeburger Allee 160
# 23562 Luebeck
# Germany
#
# http://www.imbs-luebeck.de
# wright@imbs.uni-luebeck.de
# -------------------------------------------------------------------------------

##' Grow two random forests on two cross-validation folds. 
##' Instead of out-of-bag data, the other fold is used to compute permutation importance.
##' Related to the novel permutation variable importance by Janitza et al. (2015).
##'
##' @title Hold-out random forests
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame}, \code{matrix} or \code{gwaa.data} (GenABEL).
##' @param ... Further arguments passed to ranger(). 
##' @return Hold-out random forests with variable importance.
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @references
##'   Janitza, S., Celik, E. & Boulesteix, A.-L., (2015). A computationally fast variable importance test for random forest for high dimensional data, Technical Report 185, University of Munich, \url{https://epub.ub.uni-muenchen.de/25587}. \cr
##' @export 
holdoutRF <- function(formula, data, ...) {
  ## Split data
  n <- nrow(data)
  weights <- rbinom(n, 1, 0.5)
  
  ## Grow RFs
  res <- list(
    rf1 = ranger(formula = formula, data = data, importance = "permutation",  
                 case.weights = weights, replace = FALSE, holdout = TRUE, ...),
    rf2 = ranger(formula = formula, data = data, importance = "permutation",
                 case.weights = 1-weights, replace = FALSE, holdout = TRUE, ...)
  )
  
  ## Compute importance
  res$variable.importance <- (res$rf1$variable.importance + res$rf2$variable.importance)/2
  res$treetype <- res$rf1$treetype
  res$importance.mode <- res$rf1$importance.mode
  class(res) <- "holdoutRF"
  
  res
}