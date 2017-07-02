## Test environments
* OS X 10.11.6 install, R version 3.4.0
* Ubuntu 3.13.0-100-generic, R version 3.4.0
* Windows, Microsoft R Open 3.3.2

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs of interest:

* checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘bindings’

In our build set-up, the contents of two directories (core C++, and Rcpp bindings must be copied into `src`. It was cleaner to keep the bindings in a designated directory and fully regenerate `src` on each build.

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

We need to specify a wildcard pattern for SOURCES to allow for a hierarchical structure in the `src` folder, and according to the CRAN documentation, there does not seem to be a simple cross-platform way to accomplish this: https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Compiling-in-sub_002ddirectories

## Downstream dependencies
We have no known downstream dependencies.
