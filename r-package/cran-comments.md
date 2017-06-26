## Test environments
* OS X 10.11.6 install, R version 3.4.0
* Ubuntu 3.13.0-100-generic, R version 3.4.0
* Windows, Microsoft R Open 3.3.2

## R CMD check results
There were no ERRORs. 

There was 1 WARNING:

* checking for GNU extensions in Makefiles ... WARNING
Found the following file(s) containing GNU extensions:
  src/Makevars

We need to specify a wildcard pattern for SOURCES to allow for a hierarchical structure in the `src` folder, and are not aware of a cross-platform way to accomplish this.

There was 1 NOTE of interest:

* checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘bindings’

In our build set-up, the contents of two directories (core C++, and Rcpp bindings must be copied into `src`. It was cleaner to keep the bindings in a designated directory and fully regenerate `src` on each build.

## Downstream dependencies
We have no known downstream dependencies.