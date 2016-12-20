library(Rcpp)
library(roxygen2)
library(devtools)
library(testthat)

package.name <- "gradient.forest"

compileAttributes(package.name)
roxygenize(package.name)

clean_dll(package.name)
build(package.name)
check(package.name)

install(package.name)
library(package.name, character.only = TRUE)

test_package(package.name)
