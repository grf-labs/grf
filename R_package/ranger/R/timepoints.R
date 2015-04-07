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
timepoints <- function(x, ...)  UseMethod("timepoints")

##' Extract unique death times of Ranger Survival prediction object.
##'
##'
##' @title Ranger timepoints
##' @param x Ranger Survival prediction object.
##' @param ... Further arguments passed to or from other methods.
##' @return Unique death times
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @export
timepoints.ranger.prediction <- function(x, ...) {
  if (class(x) != "ranger.predicton") {
    stop("Object ist no ranger.prediction object.")
  }
  if (x$treetype != "Survival") {
    stop("No timepoints found. Object is no Survival prediction object.")
  }
  if (is.null(x$unique.death.times)) {
    stop("No timepoints found.")
  }
  return(x$unique.death.times)
}

##' Extract unique death times of Ranger Survival forest
##'
##'
##' @title Ranger timepoints
##' @param x Ranger Survival forest object.
##' @param ... Further arguments passed to or from other methods.
##' @return Unique death times
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @aliases timepoints
##' @export
timepoints.ranger <- function(x, ...) {
  if (class(x) != "ranger") {
    stop("Object ist no ranger object.")
  }
  if (x$treetype != "Survival") {
    stop("No timepoints found. Object is no Survival forest.")
  }
  if (is.null(x$unique.death.times)) {
    stop("No timepoints found.")
  }
  return(x$unique.death.times)
}
