# Developing

In addition to providing out-of-the-box forests for quantile regression and instrumental variables, grf provides a framework for creating forests tailored to new statistical tasks. Certain components around splitting and prediction can be swapped out, within the general infrastructure for growing and predicting on trees.

## Contributing

This repository follows the standard open source protocol and setup with git where there is an abundance of existing resources to get up to speed (see for example the
contributing guidelines for well known packages in other languages, like Scikit-learn, Scipy, and pandas)

Condensed greatly, the workflow is to fork this repository, check out a branch, commit your changes (forming an ideally legible commit history),
then submitting a pull request explaining your contribution, ideally referring to the issue you created, or the issue you chose to work on.

## Working with the code

The core forest implementation is written in C++, with an R interface powered by Rcpp. We recommend using a full-powered C++ IDE such as CLion, Xcode, or Visual Studio when working with the core code.

## R package

To build the R package from source, cd into `r-package` and run `build_package.R`. Required development dependencies are listed there. This mimics the tests run when submitting a pull request. Additional online package documentation is built on Travis with [pkgdown](https://pkgdown.r-lib.org/). Usage examples in the form of R Markdown files under _grf\vignettes_ are built and rendered and the R method reference (along with the package articles) is displayed according to the layout defined in _pkgdown.yml_. To build the site locally run `pkgdown::build_site()` from the R package directory.

An alternative development workflow is to use the accompanying grf.Rproj and build and test the package with RStudio's build menu, which can be convenient
for quickly iterating C++/R code changes.

### Note for Windows users:

Symlinks in the src directory point to the core C++ and R bindings. On Windows one has to clone this repository with symlinks enabled: `git clone -c core.symlinks=true https://github.com/grf-labs/grf.git` (this command needs to be run as an administrator: right click _Command Prompt -> Run as administrator_). Caveat: the above RStudio workflow is not tested on Windows.

## Core C++

### Code structure

![GRF Architecture Diagram](https://raw.githubusercontent.com/grf-labs/grf/master/diagrams/arch_diagram.png)

The forest implementation is composed of two top-level components, [ForestTrainer](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestTrainer.h) and [ForestPredictor](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestPredictor.h).

ForestTrainer drives the tree-growing process, and has two pluggable components.
* [RelabelingStrategy](https://github.com/grf-labs/grf/blob/master/core/src/relabeling/RelabelingStrategy.h) is applied before every split, and produces a set of relabelled outcomes given the observations for a group of samples. In the case of quantile forests, for example, this strategy computes the quantiles for the group of samples, then relabels them with a factor corresponding to the quantile they belong to.
* [SplittingRule](https://github.com/grf-labs/grf/blob/master/core/src/splitting/SplittingRule.h) is called to find the best split for a particular node, given a set of outcomes. There are currently implementations for standard regression and multinomial splitting.

The trained forest produces a [Forest](https://github.com/grf-labs/grf/blob/master/core/src/forest/Forest.h) object. This can then be passed to the ForestPredictor to predict on test samples. The predictor has a pluggable 'prediction strategy', which computes a prediction given a test sample. Prediction strategies can be one of two types:
* [DefaultPredictionStrategy](https://github.com/grf-labs/grf/blob/master/core/src/prediction/DefaultPredictionStrategy.h) computes a prediction given a weighted list of training sample IDs that share a leaf with the test sample. Taking quantile forests as an example, this strategy would compute the quantiles of the weighted leaf samples.
* [OptimizedPredictionStrategy](https://github.com/grf-labs/grf/blob/master/core/src/prediction/OptimizedPredictionStrategy.h) does not predict using a list of neighboring samples and weights, but instead precomputes summary values for each leaf during training, and uses these during prediction. This type of strategy will also be passed to ForestTrainer, so it can define how the summary values are computed.

Prediction strategies can also compute variance estimates for the predictions, given a forest trained with grouped trees. Because of performance constraints, only 'optimized' prediction strategies can provide variance estimates.

A particular type of forest is created by pulling together a set of pluggable components. As an example, a quantile forest is composed of a QuantileRelabelingStrategy, ProbabilitySplittingRule, and QuantilePredictionStrategy.
The factory classes [ForestTrainers](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestTrainers.h) and [ForestPredictors](https://github.com/grf-labs/grf/blob/master/core/src/forest/ForestPredictors.h) define the common types of forests like regression, quantile, and causal forests.

### Creating a custom forest

When designing a forest for a new statistical task, we suggest following the format suggested in the paper: provide a custom RelabelingStrategy and PredictionStrategy, but use the standard RegressionSplittingRule. If modifications to the splitting rule must be made, they should be done carefully, as it is an extremely performance-sensitive piece of code.

To avoid the overhead of setting up new classes and Rcpp bindings, we provide a template for a 'custom' forest. To get started, fill in the implementation for CustomRelabelingStrategy and CustomPredictionStrategy. By default, the custom forest uses the standard regression splitting rule -- if you'd like to change this, consult the ForestTrainers.custom_trainer method.

This forest is made available in R as `custom.forest` (with `predict.custom.forest`). You can find a starter template for a test exercising the forest in `testthat/test_custom_forest.R`. Note that you'll need to re-run `build_package.R` after making changes to the C++ source.
