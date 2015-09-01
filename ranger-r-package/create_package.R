
library(Rcpp)
library(roxygen2)
library(devtools)
library(testthat)

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

## Add Documentation
roxygenize(package.name)

## Build/check/install/load package
clean_dll(package.name)
build(package.name)
##check(package.name)
install(package.name)
library(package.name, character.only = TRUE)

## Test
test_dir("test/", reporter = "summary")
