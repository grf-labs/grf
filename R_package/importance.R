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

##' Extract variable importance of Ranger object.
##'
##'
##' @title Ranger variable importance
##' @param x Ranger object.
##' @param ... Further arguments passed to or from other methods.
##' @return Variable importance measures.
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
##' @aliases importance
##' @export
importance.ranger <- function(x) {
  if (class(x) != "ranger") {
    stop("Object ist no ranger object.")
  }
  if (is.null(x$variable.importance) | length(x$variable.importance) < 1) {
    stop("No variable importance found. Please use 'importance' option when growing the forest.")
  }
  return(x$variable.importance)
}
