library(Rcpp)
library(devtools)
library(testthat)

package.name <- "grf"
package.src <- file.path(package.name, "src")
package.bindings <- file.path(package.name, "bindings")

# Copy Rcpp bindings and C++ source into the package src directory.
unlink(package.src, recursive = TRUE)
dir.create(package.src)

binding.files <- list.files(package.bindings, full.names = TRUE)
file.copy(binding.files, package.src, recursive = FALSE)

file.copy("../core/src", package.src, recursive = TRUE)
file.copy("../core/third_party", package.src, recursive = TRUE)

# Run Rcpp and build the package.
compileAttributes(package.name)
clean_dll(package.name)
document(package.name)
build(package.name)

# Test installation and run some smoke tests.
install(package.name)
library(package.name, character.only = TRUE)

test_package(package.name)
