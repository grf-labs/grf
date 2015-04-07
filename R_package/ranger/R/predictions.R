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
# Institut für Medizinische Biometrie und Statistik
# Universität zu Lübeck
# Ratzeburger Allee 160
# 23562 Lübeck
#
# http://www.imbs-luebeck.de
# wright@imbs.uni-luebeck.de
# -------------------------------------------------------------------------------

##' @export
predictions <- function(x, ...)  UseMethod("predictions")

##' Extract predictions of Ranger prediction object.
##'
##'
##' @title Ranger predictions
##' @param x Ranger prediction object.
##' @param ... Further arguments passed to or from other methods.
##' @return Predictions: Classes for Classification forests, Numerical values for Regressions forests and the estimated survival functions for all individuals for Survival forests.
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @aliases predictions
##' @export
predictions.ranger.prediction <- function(x, ...) {
  if (class(x) != "ranger.prediction") {
    stop("Object ist no ranger.prediction object.")
  }
  if (x$treetype == "Classification" | x$treetype == "Regression" | x$treetype == "Probability estimation") {
    if (is.null(x$predictions)) {
      stop("No predictions found.")
    } else {
      return(x$predictions)
    }
  } else if (x$treetype == "Survival") {
    if (is.null(x$survival)) {
      stop("No predictions found.")
    } else {
      return(x$survival)
    }
  } else {
    stop("Unknown tree type.")
  }
}

##' Extract training data predictions of Ranger object.
##'
##'
##' @title Ranger predictions
##' @param x Ranger object.
##' @param ... Further arguments passed to or from other methods.
##' @return Predictions: Classes for Classification forests, Numerical values for Regressions forests and the estimated survival functions for all individuals for Survival forests.
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
predictions.ranger<- function(x, ...) {
  if (class(x) != "ranger") {
    stop("Object ist no ranger object.")
  }
  if (x$treetype == "Classification" | x$treetype == "Regression" | x$treetype == "Probability estimation") {
    if (is.null(x$predictions)) {
      stop("No predictions found.")
    } else {
      return(x$predictions)
    }
  } else if (x$treetype == "Survival") {
    if (is.null(x$survival)) {
      stop("No predictions found.")
    } else {
      return(x$survival)
    }
  } else {
    stop("Unknown tree type.")
  }
}
