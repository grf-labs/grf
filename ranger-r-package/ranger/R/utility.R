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

##' Prediction with new data and a saved forest from Ranger.
##'
##' @title Recode factor levels in data.frame
##' @param data Data frame with factor columns to recode
##' @param response Response variable
##' @return Object of class \code{ranger.prediction} with elements
##' @seealso \code{\link{ranger}}
##' @author Marvin N. Wright
recode.factors <- function(data, response) {
  ## For each column
  lapply(data, function(x) {
    if (is.factor(x)) {
      levels <- levels(x)
    } else {
      levels <- unique(x)
    }
    
    ## Order factor levels
    means <- aggregate(response~x, FUN=mean)
    levels.ordered <- means$x[order(means$response)]
    
    ## Return reordered factor
    factor(x, levels = levels.ordered)
  })
}