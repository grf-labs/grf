##### Version 0.4.4
* Add p-values for variable importance
* Bug fixes

##### Version 0.4.3
* Add splitting by maximally selected rank statistics for survival forests
* Bug fixes

##### Version 0.4.2
* Add Windows multithreading support for new toolchain

##### Version 0.4.1
* Runtime improvement for regression forests on classification data

##### Version 0.4.0
* New CRAN version. New CRAN versions will be 0.x.0, development versions 0.x.y

##### Version 0.3.9
* Reduce memory usage of savest forest objects (changed child.nodeIDs interface)

##### Version 0.3.8
* Remove tuning functions, please use mlr or caret

##### Version 0.3.7
* Fix bug with alternative interface and prediction
* Small fixes

##### Version 0.3.6
* Add keep.inbag option to track in-bag counts
* Add option sample.fraction for fraction of sampled observations

##### Version 0.3.5
* Add tree-wise split.select.weights

##### Version 0.3.4
* Add predict.all option in predict() to get individual predictions for each tree for classification and regression
* Small changes in documentation

##### Version 0.3.3
* Add case-specific random forests

##### Version 0.3.2
* Add case weights (weighted bootstrapping or subsampling)

##### Version 0.3.1
* Catch error of outdated gcc not supporting C++11 completely

##### Version 0.3.0
* Allow the user to interrupt computation from R
* Transpose classification.table and rename to confusion.matrix
* Respect R seed for prediction
* Memory improvements for variable importance computation
* Fix bug: Probability prediction for single observations
* Fix bug: Results not identical when using alternative interface

##### Version 0.2.7 
* Small fixes for Solaris compiler

##### Version 0.2.6 
* Add C-index splitting
* Fix NA SNP handling

##### Version 0.2.5 
* Fix matrix and gwaa alternative survival interface
* Version submitted to JSS

##### Version 0.2.4 
* Small changes in documentation

##### Version 0.2.3 
* Preallocate memory for splitting

##### Version 0.2.2 
* Remove recursive splitting

##### Version 0.2.1 
* Allow matrix as input data in R version

##### Version 0.2.0 
* Fix prediction of classification forests in R

##### Version 0.1.9 
* Speedup growing for continuous covariates
* Add memory save option to save memory for very large datasets (but slower)
* Remove memory mode option from R version since no performance gain

##### Version 0.1.8 
* Fix problems when using Rcpp <0.11.4

##### Version 0.1.7 
* Add option to split on unordered categorical covariates

##### Version 0.1.6 
* Optimize memory management for very large survival forests

##### Version  0.1.5 
* Set required Rcpp version to 0.11.2
* Fix large $call objects when using BatchJobs
* Add details and example on GenABEL usage to documentation
* Minor changes to documentation

##### Version 0.1.4 
* Speedup for survival forests with continuous covariates
* R version: Generate seed from R. It is no longer necessary to set the
  seed argument in ranger calls.

##### Version 0.1.3 
* Windows support for R version (without multithreading)

##### Version 0.1.2 
* Speedup growing of regression and probability prediction forests
* Prediction forests are now handled like regression forests: MSE used for
	prediction error and permutation importance
* Fixed name conflict with randomForest package for "importance"
* Fixed a bug: prediction function is now working for probability
	prediction forests
* Slot "predictions" for probability forests now contains class probabilities
* importance function is now working even if randomForest package is
	loaded after ranger
* Fixed a bug: Split selection weights are now working as expected
* Small changes in documentation
