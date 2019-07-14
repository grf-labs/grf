# Builds and tests the GRF package.
#
# To build the package for development:
#   `Rscript build_package.R`
#
# To prepare a CRAN build:
#   `Rscript build_package.R --as-cran`

args <- commandArgs(TRUE)
library(Rcpp)
library(devtools)
library(testthat)
library(roxygen2)
library(lintr)

package.name <- "grf"
package.src <- "grf/src"

# Check code style consistency
linters <- with_defaults(
  line_length_linter(120), # Max line length = 120
  object_name_linter = NULL, # To allow variable.names + function_names
  commented_code_linter = NULL, # Misc. false positives
  object_usage_linter = NULL, # Misc. false positives
  spaces_left_parentheses_linter = NULL, # Misc. false positives
  default = default_linters);

lint_res <- lint_package(path = package.name, linters = linters, exclusions = file.path(package.name, "R/RcppExports.R"))
lint_res

if (length(lint_res) > 0) quit(status = 1)

# If built for CRAN, exlude all test except ones with "cran" in the filename
# by adding the following regex to .Rbuildignore.
if (!is.na(args[1]) && args[1] == "--as-cran") {
  write_union("grf/.Rbuildignore", "^tests/testthat/test_((?!cran).).*")
}

# Copy Rcpp bindings and C++ source into the package src directory. Note that we
# don't copy in third_party/Eigen, because for the R package build we provide
# access to the library through RcppEigen.
unlink(package.src, recursive = TRUE)
dir.create(package.src)

binding.files <- list.files("grf/bindings", full.names = TRUE)
file.copy(binding.files, package.src, recursive = FALSE)
file.copy("../core/src", package.src, recursive = TRUE)
file.copy("../core/third_party/optional", package.src, recursive = TRUE)
file.copy("../core/third_party/Distributions", package.src, recursive = TRUE)

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
