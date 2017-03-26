library(Rcpp)
library(devtools)
library(testthat)

package.name <- "gradient.forest"

compileAttributes(package.name)

clean_dll(package.name)
build(package.name)

install(package.name)
library(package.name, character.only = TRUE)

test_package(package.name)
