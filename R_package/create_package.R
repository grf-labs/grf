
library(Rcpp)
library(roxygen2)
library(devtools)
library(testthat)

setwd("~/myWork/ranger/R_package")
package.name <- "ranger"

## Set version and date
package.version <- scan("../source/src/version.h", character(0))[5]
dcf.file <- file.path(".", package.name, "DESCRIPTION")
dcf <- read.dcf(dcf.file)
dcf[1, "Version"] <- package.version
dcf[1, "Date"] <- as.character(Sys.Date())
write.dcf(dcf, dcf.file)

## Create Rcpp files
compileAttributes(package.name)

## Add globals.h to RcppExports.cpp (for building on Windows)
rcpp.exports.file <- file.path(".", package.name, "src", "RcppExports.cpp")
input <- readLines(rcpp.exports.file)
writeLines(c("#include \"globals.h\"", input), rcpp.exports.file)

## Add Documentation
roxygenize(package.name)

## Build/check/install/load package
build(package.name)
##check(package.name)
install(package.name)
library(package.name, character.only = TRUE)

## Test
test_dir("test/", reporter = "summary")

## Copy to vbox
system(paste("cp", paste(package.name, "_", package.version, ".tar.gz", sep = ""), "~/myWork/vbox_share/"))
