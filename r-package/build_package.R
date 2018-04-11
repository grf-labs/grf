library(Rcpp)
library(devtools)
library(testthat)
library(roxygen2)

package.name <- "grf"
package.src <- "grf/src"

# Copy Rcpp bindings and C++ source into the package src directory.
unlink(package.src, recursive = TRUE)
dir.create(package.src)

binding.files <- list.files("grf/bindings", full.names = TRUE)
file.copy(binding.files, package.src, recursive = FALSE)
file.copy("../core/src", package.src, recursive = TRUE)

# Auto-generate documentation files
roxygen2::roxygenise(package.name)

# Run Rcpp and build the package.
compileAttributes(package.name)
clean_dll(package.name)
build(package.name)

# Test installation and run some smoke tests.
install(package.name)
library(package.name, character.only = TRUE)
test_package(package.name)
