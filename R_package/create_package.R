
## Package infos
package.name <- "ranger"
package.title <- "A fast implementation of Random Forests"
package.description <- paste("Ranger is a fast implementation of Random Forests,",
  "particularly suited for high dimensional data.",
  "For now, ensembles of classification, regression and survival trees are supported.",
  "With Ranger data from Genome-wide association (GWA) studies can be analyzed efficiently.",
  "In addition to data frames, datasets of class 'gwaa.data' (R package GenABEL) can be directly analyzed.")
package.depends <- "R (>= 3.1)"
package.suggests <- "survival"

package.author <- "Marvin N. Wright"
package.email <- "wright@imbs.uni-luebeck.de"
package.license <- "GPL-3"

## Add C++ and R files
setwd("~/myWork/ranger/R_package")
src.file <- c("Makevars", "Makevars.win", "rangerCpp.cpp",
              list.files("../source/src", recursive = TRUE, full.names = TRUE))
src.file <- src.file[src.file != ("../source/src/main.cpp") &
                       src.file != ("../source/src/utility/ArgumentHandler.cpp") &
                       src.file != ("../source/src/utility/ArgumentHandler.h")]
code.file <- c("ranger.R", "predict.R", "print.R", "importance.R", "predictions.R", 
               "timepoints.R", "tune.R")

## Get version
package.version <- scan("../source/src/version.h", character(0))[5]

## Cleanup
remove.packages(package.name)
unlink(file.path(package.name, ""), recursive=TRUE)
unlink(paste(package.name, "_", package.version, ".tar.gz", sep = ""))
unlink(paste(package.name, ".Rcheck/", sep = ""), recursive = TRUE)

## Create package
library(Rcpp)
Rcpp.package.skeleton(package.name, attributes = TRUE, force = TRUE,
                      cpp_files = src.file, code_files = code.file,
                      example_code = FALSE, author = package.author,
                      email = package.email,
                      license = package.license)
compileAttributes(file.path(".", package.name, ""))
unlink(file.path(package.name, "man", paste(package.name, "-package.Rd", sep = "")))

## TODO: Instead use inst/include/ranger.h ?
## Add globals.h to RcppExports.cpp (for building on Windows)
# rcpp.exports.file <- file.path(".", package.name, "src", "RcppExports.cpp")
# input <- readLines(rcpp.exports.file)
# writeLines(c("#include \"globals.h\"", input), rcpp.exports.file)

## Change DESCRIPTION file
dcf.file <- file.path(".", package.name, "DESCRIPTION")
dcf <- read.dcf(dcf.file)
dcf[1, "Title"] <- package.title
dcf[1, "Description"] <- package.description
dcf[1, "Version"] <- package.version
dcf.new <- cbind(cbind(dcf, "Depends" = package.depends), "Suggests" = package.suggests)
write.dcf(dcf.new, dcf.file)

## Add Documentation
library(roxygen2)
roxygenize(package.name)

## Update NAMESPACE file
namespace.file <- file.path(".", package.name, "NAMESPACE")
cat(paste("useDynLib(", package.name, ")\n", sep = ""), file = namespace.file, append = TRUE)
cat("importFrom(Rcpp, evalCpp)\n", file = namespace.file, append = TRUE)

## Build package
system(paste("R CMD build", package.name))

## Check
##system(paste("R CMD check", package.name))

## Install and load package
install.packages(paste(package.name, "_", package.version, ".tar.gz", sep = ""), repos = NULL)
library(package.name, character.only = TRUE)

## Test
library(testthat)
test_dir("test/", reporter = "summary")
