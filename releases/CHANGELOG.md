# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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
