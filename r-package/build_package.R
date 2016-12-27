library(Rcpp)
library(devtools)
library(testthat)

package.name <- "gradient.forest"

unlink("gradient.forest/src/*")
bindings.files <- list.files("gradient.forest/bindings", recursive = TRUE, full.names = TRUE) 
core.files <- list.files("../core", "\\.(h)|(hpp)|(cpp)$", recursive = TRUE, full.names = TRUE)
file.copy(c(bindings.files, core.files), "gradient.forest/src")

compileAttributes(package.name)

clean_dll(package.name)
build(package.name)
#check(package.name)

install(package.name)
library(package.name, character.only = TRUE)

#test_package(package.name)
