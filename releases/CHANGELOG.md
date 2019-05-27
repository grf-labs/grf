# Changelog
All notable changes to grf will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [0.10.3] - 2019-05-03

### Breaking Changes
- Fix two bugs in the termination criterion for tree splitting. **IMPORTANT:** these bug fixes may cause results to change compared to previous releases, even if the same random seed is used.
  - Remove the purity condition on outcomes during splitting. For all tree types, we used to stop splitting if all outcomes in a leaf are the same. This behavior does not make sense for causal forests (which incorporates other observations besides the outcome), so it was removed. [#362](https://github.com/grf-labs/grf/pull/362).
  - Stop splitting if the objective can no longer be improved. With this change, `causal_forest` may split slightly less aggressively. [#415](https://github.com/grf-labs/grf/pull/415)

### Added
- In out-of-bag prediction, return the Monte Carlo error alongside the debiased error. [#327](https://github.com/grf-labs/grf/pull/327)
- Allow for passing a factor for the `cluster` parameter. [#329](https://github.com/grf-labs/grf/pull/329)
- Support taking a union of forests through the `merge_forests` method. [#347](https://github.com/grf-labs/grf/pull/347)
- Include a summary of the parameter tuning procedure in the forest object. [#419] (https://github.com/grf-labs/grf/pull/419)
- Add experimental support for sample weighting to regression, causal, and instrumental forests. [#376](https://github.com/grf-labs/grf/pull/376), [#418](https://github.com/grf-labs/grf/pull/418)
- Add a new forest type `boosted_regression_forest`, which applies boosting to regression forests. Allow boosting to be used during orthogonalization through the `orthog.boosting` parameter. [#388](https://github.com/grf-labs/grf/pull/388)

### Fixed
- Improve input data validation. [#354](https://github.com/grf-labs/grf/pull/354), [#378](https://github.com/grf-labs/grf/pull/378), [#430](https://github.com/grf-labs/grf/pull/430)
- Improve the `test_calibration` function by switching to one-sided p-values. [#370](https://github.com/grf-labs/grf/pull/370)
- For custom forests, fix a bug in OOB prediction where the train and tests datasets were switched. [#372](https://github.com/grf-labs/grf/pull/372)
- Decrease memory usage during training and out-of-bag prediction. [#408](https://github.com/grf-labs/grf/pull/408), [#412](https://github.com/grf-labs/grf/pull/412)
- Allow roxygen to autogenerate the `NAMESPACE` file. [#423](https://github.com/grf-labs/grf/pull/423), [#428](https://github.com/grf-labs/grf/pull/428)

## [0.10.2] - 2018-11-23
### Added
- Add support for confidence intervals in local linear regression forests.

### Changed
- Allow samples_per_cluster to be larger than smallest cluster size.

### Fixed
- Make sure average effect estimation doesn't error on data with a single feature.
- Fix a bug in local linear prediction where the penalty wasn't properly calculated.
- Fix two issues in causal forest tuning that could lead to unstable results.
- Ensure that the ATE and APE functions correctly account for cluster membership.

## [0.10.1] - 2018-09-23
### Added
- Add basic support for tree plotting (through the `plot` method).
- Add the method `test_calibration`, which performs an omnibus test for presence of heterogeneity via calibration.
- For local linear regression forests, add support for selecting the value of `ll.lambda` through cross-validation.
- Introduce a training option `honesty.fraction` that can be used to specify the fraction of data that should be used in selecting splits vs. performing estimation. Note that this parameter is only relevant when honesty is enabled (the default).
- Start a practical guide to the GRF algorithm (https://github.com/grf-labs/grf/blob/master/REFERENCE.md).
- In `average_treatment_effect` and `average_partial_effect`, add an option `subset` to support estimating the treatment effect over a subsample of the data.

### Fixed
- Fix a bug in random sampling where features listed earlier in the data matrix were more likely to be selected for splitting.
- Make sure that the sample indices returned in `get_tree` are 1-indexed, as opposed to 0-indexed.

## [0.10.0] - 2018-05-08
### Added
- Replace the local.linear option with linear.correction.variables, which allows for a subset of
  variables to be considered during local linear regression.
- Update causal_forest interface to allow user to specify Y.hat and W.hat. Note that these options supercede
  the precompute.nuisance parameter, which has been removed. To recreate the behavior of
  precompute.nuisance = TRUE, NULL can be provided for Y.hat and W.hat, and for precompute.nuisance = FALSE,
  Y.hat and W.hat should be 0.
- Add a causal forest example with variable selection and parameter tuning.

### Changed
- Adjust the defaults for the causal forest tuning algorithm.

### Fixed
- Prevent tuning down to a min.node.size of 0.

## [0.9.6] - 2018-04-13

### Added
- Debiased error criterion for measuring the out-of-bag accuracy of a forest using only a few trees.
- Automated tuning via cross-validation for regression and causal forests.
- Estimation of average partial effects with a continuous treatment.
- Overlap-weighted average treatment effects.
- Cluster-robust standard errors for regression and causal forests, and average effect estimates (contributed by @lminer).
- Locally linear prediction in regression forests (contributed by @rinafriedberg).
- Regularize splits in causal/instrumental forests via a variance penalty.

### Changed
- Avoid causal forest leaves with all treated or all control samples (controlled via stabilize.splits = TRUE).
- Store in-bag rather than out-of-bag samples to save memory.
- Only support sampling with replacement (as some features are ambiguously defined with bootstrapping).

### Fixed
- Bug in IV CI construction (#209).

## [0.9.5] - 2018-01-04

### Added
- Create a simple method for variable importance based on split frequency and depth.
- Add support for sparse data matrices of type 'dgCMatrix'.

### Changed
- Use RcppEigen for the R package, as opposed to the eigen source.

### Fixed
- Fix a few places where we were still using the old default for mtry. This issue was
causing poor performance for even moderately large numbers of features.

## [0.9.4] - 2017-11-25

### Changed
- Update the default for mtry to sqrt(p) + 20.

### Fixed
- Fix an issue where split_frequencies fails when p = 1.
- Use a Solaris-compatible version of std::sqrt.

## [0.9.3] - 2017-07-20

### Fixed
- Fix an out of bounds error when there are fewer trees than threads.
- Fix a bug in the get_tree function where the same tree was always returned.

## [0.9.2] - 2017-07-06

### Added
- Add an experimental regularized version of the regression splitting rule.

### Fixed
- Several bugfixes for CRAN compatibility.

## [0.9.1] - 2017-06-27

### Added
- Add a regression.splits option to quantile forests to allow emulating the approach in Meinhausen (2006).

## [0.9.0] - 2017-06-25
First official (beta) release. The package currently supports
- standard regression forests
- causal, instrumental, and quantile forests
- confidence intervals for causal, instrumental, and regression forests
- training 'honest' versions of the above forests
