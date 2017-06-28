# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [0.9.1] - 2017-06-27

### Added
- Add a regression.splits option to quantile forests to allow emulating the approach in Meinhause (2006).

## [0.9.0] - 2017-06-25
First official (beta) release. The package currently supports
- standard regression forests
- causal, instrumental, and quantile forests
- confidence intervals for causal, instrumental, and regression forests
- training 'honest' versions of the above forests
