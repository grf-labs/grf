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

##' Print contents of Ranger object.
##'
##'
##' @title Print Ranger
##' @param x Object of class 'ranger'.
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
print.ranger <- function(x, ...) {
  cat("Ranger result\n\n")
  cat("Call:\n", deparse(x$call), "\n\n")
  cat("Type:                            ", x$treetype, "\n")
  cat("Number of trees:                 ", x$num.trees, "\n")
  cat("Sample size:                     ", x$num.samples, "\n")
  cat("Number of independent variables: ", x$num.independent.variables, "\n")
  cat("Mtry:                            ", x$mtry, "\n")
  cat("Target node size:                ", x$min.node.size, "\n")
  cat("Variable importance mode:        ", x$importance.mode, "\n")
  if (x$treetype == "Survival") {
    cat("Number of unique death times:    ", length(x$unique.death.times), "\n")
  }
  if (x$treetype == "Classification") {
    cat("OOB prediction error:            ", sprintf("%1.2f %%", 100*x$prediction.error), "\n")
  } else {
    cat("OOB prediction error:            ", x$prediction.error, "\n")
  }
  if (x$treetype == "Regression") {
    cat("R squared:                       ", x$r.squared, "\n")
  }
}

##' Print contents of Ranger forest object.
##'
##'
##' @title Print Ranger forest
##' @param x Object of class 'ranger.forest'.
##' @param ... further arguments passed to or from other methods.
##' @author Marvin N. Wright
##' @export
print.ranger.forest <- function(x, ...) {
  cat("Ranger forest object\n\n")
  cat("Type:                         ", x$treetype, "\n")
  cat("Number of trees:              ", x$num.trees, "\n")
  if (x$treetype == "Survival") {
    cat("Number of unique death times: ", length(x$unique.death.times), "\n")
  }
}

##' Print contents of Ranger prediction object.
##'
##'
##' @title Print Ranger prediction
##' @param x Object of class 'ranger.prediction'.
##' @param ... further arguments passed to or from other methods.
##' @author Marvin N. Wright
##' @export
print.ranger.prediction <- function(x, ...) {
  cat("Ranger prediction\n\n")
  cat("Type:                            ", x$treetype, "\n")
  cat("Sample size:                     ", x$num.samples, "\n")
  cat("Number of independent variables: ", x$num.independent.variables, "\n")
  if (x$treetype == "Survival") {
    cat("Number of unique death times:    ", length(x$unique.death.times), "\n")
  }
}

str.ranger.forest <- function(object, max.level = 2, ...) {
  class(object) <- "list"
  str(object, max.level = max.level, ...)
}

str.ranger <- function(object, max.level = 2, ...) {
  class(object) <- "list"
  str(object, max.level = max.level, ...)
}
